import os
import pandas as pd
import numpy as np
from functools import lru_cache
import dalmatian
from scipy.stats import beta

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

    final_fig = make_subplots(rows=1,cols=2, shared_yaxes=True, column_widths=[0.77, 0.25])
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

    return mod_tab_df, maf_file, mut_annot_list

def calculate_multiplicity(maf_df, alpha):
    # Local copy number
    q_values = maf_df['q_hat']

    # calculates Allele fraction
    allele_fraction = maf_df['alt'] / (maf_df['alt']+maf_df['ref'])
    multiplicities = allele_fraction * (alpha*q_values + (1-alpha)*2) / alpha

    return(multiplicities)

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

def add_precalled_purities_to_pairs(pairs_df, sample_df, precalled_purity_col_nm='PCA__ABSOLUTE__Cancer_DNA_fraction'):
    """
    Makes a new dataframe with the precalled purity values from the sample_df added to the columns from the pairs_df

    Parameters
    ==========
    pairs_df: pd.DataFrame, 
        dataframe of the pairs sample data

    sample_df: pd.DataFrame, 
        dataframe of the sample data
    
    precalled_purity_col_nm: str,
        name of the column with the precalled purity values in the sample_df dataframe

    Returns
    =======
        Merged dataframe with the precalled purity values added as a new column to the original pairs_df dataframe
    """

    pairs_df = pairs_df.copy()
    sample_df = sample_df.copy()

    precalled_purities_df = sample_df[[precalled_purity_col_nm]]
    precalled_purities_df = precalled_purities_df[precalled_purities_df[precalled_purity_col_nm].notna()]

    pairs_with_precalled_purity_df = pairs_df.merge(precalled_purities_df, how="inner", left_on="case_sample", right_on="sample_id")
    pairs_with_precalled_purity_df = pairs_with_precalled_purity_df.rename(columns={precalled_purity_col_nm:'precalled_purity_values'})

    return pairs_with_precalled_purity_df

def gen_mult_allele_subplots(
    multiplicity_fig,
    allele_fraction_fig,
    csize=None
):
    """
    Generates the plotly subplot for the multiplicity and allele fraction graphs

    Parameters
    ==========
    multiplicity_fig: go.Figure
        figure with the multiplicity data

    allele_fraction_fig: go.Figure
        figure with the allele fraction data

    Return
    ======
    go.Figure    
        subplot figure with both the multiplicity (col 1) and allele fraction (col 2) traces on it 
    """
    if csize is None:
        csize = CSIZE_DEFAULT

    multiplicity_allele_subplots_fig = go.Figure()

    multiplicity_allele_subplots_fig = make_subplots(rows=1, cols=2, column_widths=[0.79, 0.25],
                                                     subplot_titles=('Multiplicity Plot', 'Mutation Allele Fraction Plot'))

    # Add traces from the first figure, mut_fig_with_lines, to the left subplot (row=1, col=1)
    for trace in multiplicity_fig.data:
        multiplicity_allele_subplots_fig.add_trace(trace, row=1, col=1)

    # Add traces from the second figure, allele_fraction_fig, to the right subplot (row=1, col=2)
    for trace in allele_fraction_fig.data:
        multiplicity_allele_subplots_fig.add_trace(trace, row=1, col=2)

    add_background(multiplicity_allele_subplots_fig, csize.keys(), csize, height=100, plotly_row=1, plotly_col=1)
    multiplicity_allele_subplots_fig.update_xaxes(multiplicity_fig.layout.xaxis, row=1, col=1)
    multiplicity_allele_subplots_fig.update_yaxes(multiplicity_fig.layout.yaxis, row=1, col=1)

    for yval in [1,2]:
        multiplicity_allele_subplots_fig.add_hline(y=yval, row=1,col=1, line_dash="dash", line_color='black', line_width=1)
            
    multiplicity_allele_subplots_fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", showlegend=False)
    multiplicity_allele_subplots_fig.update_xaxes(allele_fraction_fig.layout.xaxis, row=1, col=2)
    multiplicity_allele_subplots_fig.update_yaxes(allele_fraction_fig.layout.yaxis, row=1, col=2)

    return multiplicity_allele_subplots_fig

def gen_mut_allele_fraction_plot(
        maf_df, 
    ):
    """ 
    Generates the mutation allele fraction plot based on the maf file from absolute solutions report

    Parameters
    ==========
    maf_df: pd.DataFrame
        dataframe with the maf file mutation data

    Return
    ======
        tuple
            plotly.Figure
                plot for the allele fraction for each mutation
    """
    clonal_probabilities = maf_df['Pr_somatic_clonal'].values
    ssnv_skew = maf_df['SSNV_skew'].values[0]
    alpha = maf_df['purity'].values[0]
    line_colors = ["blue", "grey"]
    n_grid = 300
    allele_fraction_beta_distributions = get_mut_beta_densities(maf_df, n_grid)

    # Creates a normalized distribition for plotting
    normalized_values = np.linspace(1, n_grid, n_grid) / (n_grid + 1)
    normalized_values_matrix = np.tile(normalized_values, (allele_fraction_beta_distributions.shape[0], 1))  # Repeat normalized distribution for each mutation
    
    # creates a figure with the probability distribution (individual and total subclonal and clonal)
    fig = go.Figure()
    fig = draw_mut_beta_densities(allele_fraction_beta_distributions, clonal_probabilities, line_colors)

    # Add vertical lines for alpha / 2 and alpha * SSNV_skew / 2
    fig.add_shape(
        type='line',
        x0=alpha / 2,
        x1=alpha / 2,
        y0=0,
        y1=1,
        line=dict(color="green", dash='dot', width=2),
        name="alpha / 2"
    )
    fig.add_annotation(
        x=alpha / 2,
        y=1.15,
        text='alpha / 2',
        showarrow=False,
        font=dict(color="black")
    )

    fig.add_shape(
        type='line',
        x0=alpha * ssnv_skew / 2,
        x1=alpha * ssnv_skew / 2,
        y0=0,
        y1=1,
        line=dict(color='grey', dash='dot', width=2),
        name="alpha * f_s / 2"
    )
    fig.add_annotation(
        x=alpha * ssnv_skew / 2,
        y=1.05,
        text='f_s * alpha / 2',
        showarrow=False,
        font=dict(color='black')
    )

    # Update layout of the figure
    fig.update_layout(
        title="Mutation Allele Fraction Plot",
        yaxis=dict(title="Density"),
        template="plotly_white",
        xaxis_range=[0, 1],
        yaxis_range=[0, 1.2],
        xaxis=dict(
            title="Fraction of Alternate Reads",
            tickvals=[i/10 for i in range(0, 11, 2)],  # Set the tick positions on the x axis
            ticktext=[f"{i/10}" for i in range(0, 11, 2)]),  # Set the labels for those tick positions on the x axis
        showlegend=False,
        width=400,
    )

    return fig

def draw_mut_beta_densities(
        af_beta_distributions, 
        clonal_probabilities, 
        line_colors
    ):
    """ 
    Creates the mutation beta densities lines on the plotly figure

    Parameters
    ==========
    af_beta_distributions: 2D numpy array
        beta distribution for each mutations allele fraction 

    clonal_probabilities: 1D numpy array
        probability of a mutation being clonal

    line_colors: list ["blue", "grey]
        colors for the total clonal and subclonal lines

    Return
    ======
        plotly.Figure
            plot for the allele fraction for each mutation
    """
    n_grid = af_beta_distributions.shape[1]
    normalized_values = np.arange(1, n_grid + 1) / (n_grid + 1)  # normalizing based on the number of columns in the beta distributions
    
    subclonal_probabilities = 1 - clonal_probabilities
    subclonal_probabilities[subclonal_probabilities < 0] = 0  # Fix round-off error
    
    # Set up the plot
    fig = go.Figure()

    clonal_probabilities_vector = np.reshape(clonal_probabilities, (clonal_probabilities.shape[0], 1))
    subclonal_probabilities_vector = np.reshape(subclonal_probabilities, (subclonal_probabilities.shape[0], 1))

    clonal_probabilities_vector = np.nan_to_num(clonal_probabilities_vector)
    subclonal_probabilities_vector = np.nan_to_num(subclonal_probabilities_vector)
    af_beta_distributions = np.nan_to_num(af_beta_distributions)

    # calculates the clonal and subclonal distributions
    clonal_distributions = af_beta_distributions * clonal_probabilities_vector
    subclonal_distributions = af_beta_distributions * subclonal_probabilities_vector
    
    # Plots the normalized total sum of subclonal distributions
    fig.add_trace(go.Scatter(
        x=normalized_values,
        y=np.sum(subclonal_distributions, axis=0) / np.max(np.sum(subclonal_distributions, axis=0)), # normalized sum of each subclonal distribution 
        mode='lines',
        line=dict(color=line_colors[0], dash='dash'),
        name="Subclonal Total"
    ))

    # Plots the normalized total sum of clonal distributions
    fig.add_trace(go.Scatter(
        x=normalized_values,
        y=np.sum(clonal_distributions, axis=0) / np.max(np.sum(clonal_distributions, axis=0)), # normalized sum of each clonal distribution 
        mode='lines',
        line=dict(color=line_colors[1], dash='dash'),
        name="Clonal Total"
    ))
    
    return fig

def get_mut_beta_densities(
        maf_df, 
        n_cols=100
    ):
    """
    Calculates the beta distributions for each mutation

    Parameters
    ==========
    maf_df: pd.DataFrame
        maf mutation data

    n_cols: int
        number of columns in the beta distribution matrix 

    Return
    ======
    2D numpy array
    """
    # Initialize mutation grid with zeros
    mut_grid = np.zeros((maf_df.shape[0], n_cols))
    
    # Calculate coverage and allele frequency
    total_num_reads = maf_df["alt"] + maf_df["ref"]
    allele_fraction = maf_df["alt"] / total_num_reads
    uniform_vals = np.linspace(0, 1, n_cols)
    
    # Calculate beta densities for each mutation
    for i in range(maf_df.shape[0]):
        alpha_parameter = total_num_reads[i] * allele_fraction[i] + 1 
        beta_parameter = total_num_reads[i] * (1 - allele_fraction[i]) + 1
        mut_grid[i, :] = beta.pdf(uniform_vals, alpha_parameter, beta_parameter) * 1 / n_cols  # calculates and normalizes the beta distribution 
    
    # Check for NaN values
    if np.any(np.isnan(mut_grid)):
        raise ValueError("NaN values detected in the mutation grid")
    
    return mut_grid