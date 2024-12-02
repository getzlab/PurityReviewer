import os
import pandas as pd
import numpy as np
from functools import lru_cache
import dalmatian

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects
from cnv_suite.visualize import plot_acr_interactive, add_background, update_cnv_scatter_sigma_toggle, prepare_df
from cnv_suite import calc_avg_cn
from natsort import natsorted

from AnnoMate.AppComponents.utils import freezeargs, cached_read_csv


CSIZE_DEFAULT = {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260,
                 '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747,
                 '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
                 '16': 90354753, '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
                 '21': 48129895, '22': 51304566, '23': 156040895, '24': 57227415}


def get_cum_sum_csize(csize):
    """
    Gets the cumulative positions of chromosomes from a dictionary of chromosome lengths

    Parameters
    ==========
    csize: Dict
        Ordered dictionary with chromosome names as the key and the chromosome's length as the value

    Returns
    =======
    Dict
        Dictionary with chromosome names as the key and the chromosome's cumulative end position as the value
    """
    cum_sum_csize = {}
    cum_sum = 0
    for chrom, size in csize.items():
        cum_sum_csize[chrom] = cum_sum
        cum_sum += size
    return cum_sum_csize


def plot_cnp_histogram(
    seg_df,
    mu_major_col,
    mu_minor_col,
    length_col,
    max_mu=2.0,
    step=0.05,
    fig=None,
    fig_row=None,
    fig_col=None,
):
    """
    Generates marginal histogram of segments over allelic copy ratio bins
    
    Parameters
    ==========
    seg_df: pd.DataFrame
        
    mu_major_col: str
        Column in seg_df corresponding to the mean allelic copy ratio of the major allele
        
    mu_minor_col: str
        Column in seg_df corresponding to the mean allelic copy ratio of the minor allele
        
    length_col: str
        Column in seg_df corresponding to the length of the segment
        
    max_mu: float, default=2.0
        maximum mu value to plot
        
    step: float, default=0.05
        bin size for mu values to calculate the histogram
        
    fig: plotly.Figure
        Plotly figure
        
    fig_row: int
        1-indexed row index for the corresponding subplot in fig
        
    fig_col: int
        1-indexed column index for the corresponding subplot in fig
        
    Returns
    =======
    plotly.Figure
        Histogram of frequency of genomic material across allelic copy ratio bins.
    """
    # bin over 
    mu_bins = np.arange(0, max_mu + step, step)
    mu_bins_counts = {b: 0 for b in mu_bins}
    for _, r in seg_df.iterrows():
        mu_bins_counts[mu_bins[np.argmax(r[mu_major_col] < mu_bins) - 1]] += r[length_col]
        mu_bins_counts[mu_bins[np.argmax(r[mu_minor_col] < mu_bins) - 1]] += r[length_col]
    
    # add half step to make them centered around the bin
    half_step = step / 2.0
    mu_bins_counts = {mu_bin + half_step: val for mu_bin, val in mu_bins_counts.items()}
    
    mu_bin_counts_df = pd.DataFrame.from_dict(mu_bins_counts, 
                                              orient='index',
                                              columns=['count'])
    mu_bin_counts_df.index.name = 'mu_bin'
    mu_bin_counts_df = mu_bin_counts_df.reset_index()
    if fig is None:
        fig = go.Figure()
    bar_trace = go.Bar(x=mu_bin_counts_df['count'], y=mu_bin_counts_df['mu_bin'], orientation='h')
    if fig_row is not None:
        fig.add_trace(bar_trace, row=fig_row, col=fig_col)
        fig.update_xaxes(title_text='Length Count',
                         row=fig_row, col=fig_col)
    else:
        fig.add_trace(bar_trace)
        fig.update_xaxes(title_text='Length Count')

    fig.update_layout(showlegend=False)
    return fig


#@freezeargs
#@lru_cache(maxsize=32)
def gen_mut_figure(maf_df,
                   chromosome_col='Chromosome',
                   start_position_col='Start_position',
                   hugo_symbol_col='Hugo_Symbol',
                   variant_type_col='Variant_Type',
                   alt_count_col='t_alt_count',
                   ref_count_col='t_ref_count',
                   hover_data=None,
                   csize=None
                  ):
    """
    Generates plot with genome on the X axis and VAF on the y axis, and plots mutations from a maf file. 

    Parameters
    ==========
    maf_fn: str
        Column in the data df with each sample's path to its maf file. This maf file should be filtered to only have somatic mutations.
        
    chromosome_col: str, default="Chromosome"
        Column in the maf file corresponding to the chromosome the mutation is in
        
    start_position_col: str, default="Start_position"
        Column in the maf file corresponding to the start position of the mutation
        
    hugo_symbol_col: str, default="Hugo_Symbol"
        Column in the maf file corresponding to the name of the gene that a mutation is in
        
    variant_type_col: str, default="Variant_Type"
        Column in the maf file corresponding to the variant type of the mutation (SNP/INDEL)
        
    alt_count_col: str, default="t_alt_count"
        Column in the maf file corresponding to the alt count of the mutation
        
    ref_count_col: str, default="t_ref_count"
        Column in the maf file corresponding to the ref count of the mutation
        
    hover_data: List[str]
        Additional columns in the maf file to display when hovering over a mutation in the plot
        
    csize: Dict
        Dictionary with chromosome sizes (see AppComponents.utils.CSIZE_DEFAULT for hg19)

    Returns
    =======
    plotly.Figure object 
    """

    hover_data = [] if hover_data is None else list(hover_data)
    csize = CSIZE_DEFAULT if csize is None else csize
    cum_sum_csize = get_cum_sum_csize(csize)

    chr_order = natsorted(list(csize.keys()))
    chrom_start = {chrom: start for (chrom, start) in
                   zip(np.append(chr_order, 'Z'), np.cumsum([0] + [csize[a] for a in chr_order]))}

    #maf_df = cached_read_csv(maf_fn, sep='\t', encoding='iso-8859-1')
    if maf_df[chromosome_col].dtype == 'object':
        maf_df[chromosome_col].replace({'X': 23, 'Y': 24}, inplace=True)
    maf_df[chromosome_col] = maf_df[chromosome_col].astype(str)
    
    maf_df['new_position'] = maf_df.apply(lambda r: cum_sum_csize[r[chromosome_col]] + r[start_position_col], axis=1)

    #maf_df['tumor_f'] = maf_df[alt_count_col] / (maf_df[alt_count_col] + maf_df[ref_count_col])
    
    fig = px.scatter(maf_df, x='new_position', y='multiplicity', marginal_y='histogram', hover_data=hover_data)

    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')

    fig.update_xaxes(
        showgrid=False,
        zeroline=False,
        tickvals=np.asarray(list(chrom_start.values())[:-1]) + np.asarray(list(csize.values())) / 2,
        ticktext=chr_order,
        tickfont_size=10,
        tickangle=0,
        range=[0, chrom_start['Z']]
    )
    fig.update_xaxes(title_text="Chromosome")
                      
                      
    fig.update_yaxes(range=[0, 2.5])

    final_fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.77, 0.25])
    for t in fig.data:
        final_fig.add_trace(t, row=1, col=1)

    add_background(final_fig, csize.keys(), csize, height=100, plotly_row=1, plotly_col=1)
    final_fig.update_xaxes(fig.layout.xaxis, row=1, col=1)
    final_fig.update_yaxes(fig.layout.yaxis, row=1, col=1)
    final_fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", showlegend=False)
    
    return final_fig


def gen_cnp_figure(acs_fn,
                   sigmas=True, 
                   mu_major_col='mu.major', 
                   mu_minor_col='mu.minor', 
                   length_col='length',
                   csize=None
                  ):
    """
    Parameters
    ==========
    acs_fn: str, Path
        path to AllelicCapSeg-like CNV file
        
    sigmas: bool
        To plot the sigmas or not
        
    mu_major_col: str
        column in the acs_fn table corresponding to the allelic copy ratio of the major allele (mu.major)
        
    mu_minor_col: str
        column in the acs_fn table corresponding to the allelic copy ratio of the minor allele (mu.minor)
        
    length_col: str
        column in the acs_fn table corresponding to the length of the segment
        
    csize: Dict
        Dictionary with chromosome sizes (see AppComponents.utils.CSIZE_DEFAULT for hg19)
    
    Returns
    =======
    plotly.Figure
        Plot of the alellic copy ratios across the genome
    """
    csize = CSIZE_DEFAULT if csize is None else csize
    cnp_fig = _gen_cnp_figure_cache(acs_fn, mu_major_col, mu_minor_col, length_col, csize)
    if not sigmas:
        update_cnv_scatter_sigma_toggle(cnp_fig, sigmas)

    return cnp_fig


@freezeargs
@lru_cache(maxsize=32)
def _gen_cnp_figure_cache(acs_fn,
                          mu_major_col,
                          mu_minor_col,
                          length_col,
                          csize):
    """
    Cached version of the plotting function to generate the copy number profile plot and a marginal allelic copy ratio histogram

    Parameters
    ==========
    acs_fn: str, Path
        path to AllelicCapSeg-like CNV file
        
    sigmas: bool
        To plot the sigmas or not
        
    mu_major_col: str
        column in the acs_fn table corresponding to the allelic copy ratio of the major allele (mu.major)
        
    mu_minor_col: str
        column in the acs_fn table corresponding to the allelic copy ratio of the minor allele (mu.minor)
        
    length_col: str
        column in the acs_fn table corresponding to the length of the segment
        
    csize: Dict
        Dictionary with chromosome sizes (see AppComponents.utils.CSIZE_DEFAULT for hg19)
        
    Returns
    =======
    plotly.Figure
        Interactive copy number profile with marginal histogram
    """
    seg_df = cached_read_csv(acs_fn, sep='\t', encoding='iso-8859-1')

    # normalize seg file (important to do before estimating ploidy)
    seg_df[['tau','sigma.tau','mu.minor','sigma.minor','mu.major','sigma.major']] = seg_df[['tau','sigma.tau','mu.minor','sigma.minor','mu.major','sigma.major']] / calc_avg_cn(seg_df)

    acr_fig, _, _, _ = plot_acr_interactive(seg_df, csize)

    hist_fig = plot_cnp_histogram(
        seg_df,
        mu_major_col=mu_major_col,
        mu_minor_col=mu_minor_col,
        length_col=length_col
    )

    cnp_fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.77, 0.25])
    for t in acr_fig.data:
        cnp_fig.add_trace(t, row=1, col=1)

    add_background(cnp_fig, csize.keys(), csize, height=100, plotly_row=1, plotly_col=1)
    cnp_fig.update_xaxes(acr_fig.layout.xaxis, row=1, col=1)
    cnp_fig.update_yaxes(acr_fig.layout.yaxis, row=1, col=1)

    for t in hist_fig.data:
        cnp_fig.add_trace(t, row=1, col=2)

    cnp_fig.update_xaxes(hist_fig.layout.xaxis, row=1, col=2)
    cnp_fig.update_yaxes(hist_fig.layout.yaxis, row=1, col=2)
    cnp_fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", showlegend=False)

    return cnp_fig

def parse_absolute_soln_simulatedTumorData(rdata_path: str) -> pd.DataFrame:
    """
    function to convert an rdata file from ABSOLUTE into a pandas dataframe
    
    Parameters
    ==========
    rdata_path: str, Path
        LOCAL path to the rdata. This cannot read a gsurl or other remote path
    
    Returns
    =======
    pd.DataFrame
        Pandas dataframe of ABSOLUTE purity/ploidy solutions
    """
    
    absolute_rdata_cols = ['alpha', 'tau', 'tau_hat', '0_line', '1_line',
                           'sigma_H', 
                           'theta_Q', 
                           'lambda',  
                           'SCNA_likelihood', 
                           'Kar_likelihood', 
                           'SSNVs_likelihood']
    pandas2ri.activate()
    r_list_vector = robjects.r['load'](rdata_path)
    r_list_vector = robjects.r[r_list_vector[0]]

    mode_res_idx = r_list_vector.names.index('mode.res')
    sample_name_idx = r_list_vector.names.index('sample.name')
    mode_res = r_list_vector[mode_res_idx]
    mode_tab_idx = mode_res.names.index('mode.tab')
    mode_tab = mode_res[mode_tab_idx]
    mod_tab_df = pd.DataFrame(columns=absolute_rdata_cols)
    mod_tab_df['alpha'] = mode_tab[:, 0]
    mod_tab_df['tau'] = mode_tab[:, 1]
    mod_tab_df['tau_hat'] = mode_tab[:, 7]
    mod_tab_df['0_line'] = mode_tab[:, 3]
    mod_tab_df['step_size'] = mode_tab[:, 4] * 2
    mod_tab_df['1_line'] = mod_tab_df['step_size'] + mod_tab_df['0_line']
    mod_tab_df['sigma_H'] = mode_tab[:, 8]
    mod_tab_df['theta_Q'] = mode_tab[:, 11]
    mod_tab_df['lambda'] = mode_tab[:, 12]
    mod_tab_df['SCNA_likelihood'] = mode_tab[:, 15]
    mod_tab_df['Kar_likelihood'] = mode_tab[:, 17]
    mod_tab_df['SSNVs_likelihood'] = mode_tab[:, 20]
    return mod_tab_df

def parse_absolute_soln(rdata_path: str) -> pd.DataFrame: # has to be a local path   
    """
    function to convert an rdata file from ABSOLUTE into a pandas dataframe
    
    Parameters
    ==========
    rdata_path: str, Path
        LOCAL path to the rdata. This cannot read a gsurl or other remote path
    
    Returns
    =======
    pd.DataFrame
        Pandas dataframe of ABSOLUTE purity/ploidy solutions
    """
    
    absolute_rdata_cols = ['alpha', 'tau', 'tau_hat', '0_line', '1_line',
                       'sigma_H', 
                       'theta_Q', 
                       'lambda',  
                       'SCNA_likelihood', 
                       'Kar_likelihood', 
                       'SSNVs_likelihood']
    pandas2ri.activate()
    r_list_vector = robjects.r['load'](rdata_path)
    r_list_vector = robjects.r[r_list_vector[0]]
    r_data_id = r_list_vector.names[0]

    rdata_tables = r_list_vector.rx2(str(r_data_id))
    
    mode_res = rdata_tables.rx2('mode.res')
    mode_tab = mode_res.rx2('mode.tab')
    mod_tab_df = pd.DataFrame(columns=absolute_rdata_cols)
    mod_tab_df['alpha'] = mode_tab[:, 0]
    mod_tab_df['tau'] = mode_tab[:, 1]
    mod_tab_df['tau_hat'] = mode_tab[:, 7]
    mod_tab_df['0_line'] = mode_tab[:, 3]
    mod_tab_df['step_size'] = mode_tab[:, 4] * 2
    mod_tab_df['1_line'] = mod_tab_df['step_size'] + mod_tab_df['0_line']
    mod_tab_df['sigma_H'] = mode_tab[:, 8]
    mod_tab_df['theta_Q'] = mode_tab[:, 11]
    mod_tab_df['lambda'] = mode_tab[:, 12]
    mod_tab_df['SCNA_likelihood'] = mode_tab[:, 15]
    mod_tab_df['Kar_likelihood'] = mode_tab[:, 17]
    mod_tab_df['SSNVs_likelihood'] = mode_tab[:, 20]

    # Maf file
    maf_file = rdata_tables.rx2('mut.cn.dat')
    mut_annot_list = mode_res.rx2('modeled.muts')

    return mod_tab_df,maf_file,mut_annot_list

def calculate_multiplicity(maf,alpha):
    # Local copy number
    q = maf['q_hat']
    # Allele fraction
    af = maf['alt']/(maf['alt']+maf['ref'])

    # af = (alpha * mult) / (alpha * q + (1-alpha)*2)
    mult = af * (alpha*q + (1-alpha)*2) / alpha

    return(mult)

def validate_purity(x):
    return (x >=0) and (x <= 1)


def validate_ploidy(x):
    return (x >=0)


def download_data(file_to_download_path, full_local_path):
    """
    Downloads data from remote url

    Parameters
    ==========
    file_to_download_path: str
        remote url to file to download
        
    full_local_path: str
        local path to download data
        
    """
    dalmatian.getblob(file_to_download_path).download_to_filename(full_local_path)


def download_rdata(rdata_fn_s, rdata_dir, force_download=False):
    """
    Download ABSOLUTE rdata locally
    
    Parameters
    ==========
    rdata_fn_s: pd.Series
        Pandas series with sample/pair ids in the index and lists their corresponding rdata remote urls
        
    rdata_dir: str, Path
        Local path to directory to save local rdata files
        
    force_download: bool
        To force redownload even if the file already exists locally (by file name)

    Returns
    =======
    pd.Series
        Pandas series with sample/pair ids in the index and lists their corresponding rdata local urls
    """
    if not os.path.isdir(rdata_dir):
        os.mkdir(rdata_dir)

    local_rdata_dict = {}
    for pair_id, rdata_fn in rdata_fn_s.items():
        absolute_rdata_gsurl = rdata_fn

        local_absolute_rdata_fn = f'{rdata_dir}/{pair_id}.absolute.rdata'
        if not os.path.exists(local_absolute_rdata_fn) or force_download:
            download_data(absolute_rdata_gsurl, local_absolute_rdata_fn)

        local_rdata_dict[pair_id] = local_absolute_rdata_fn

    return pd.Series(local_rdata_dict)
    