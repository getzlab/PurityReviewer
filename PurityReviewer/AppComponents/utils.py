import os
import pandas as pd
import numpy as np
from functools import lru_cache
import dalmatian
from scipy.stats import beta
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects
from cnv_suite.visualize import plot_acr_interactive, add_background, update_cnv_scatter_sigma_toggle, prepare_df
from cnv_suite import calc_avg_cn
from natsort import natsorted
from scipy.interpolate import interp1d

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

def calculate_multiplicity(maf_df,alpha):
    # Local copy number
    q = maf_df['q_hat']

    # calculates Allele fraction
    allele_fraction = maf_df['alt']/(maf_df['alt']+maf_df['ref'])

    # af = (alpha * mult) / (alpha * q + (1-alpha)*2)
    mult = allele_fraction * (alpha*q + (1-alpha)*2) / alpha

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

def gen_allele_fraction_figure(maf_df):
    print("Inside the gen_allele_fraction_figure function")
    clonal_maf_df = maf_df[ maf_df['Pr_somatic_clonal'] == 'clonal' ]
    subclonal_maf_df = maf_df[ maf_df['Pr_somatic_clonal'] != 'clonal' ]

    clonal_x_values = np.linspace(0, 1, len(clonal_maf_df))
    clonal_y_values = np.array([])
    subclonal_x_values = np.linspace(0, 1, len(subclonal_maf_df))
    subclonal_y_values = np.array([])

    fig = go.Figure()

    return fig

# START OF THE PYTHON CODE VERSION OF THE ABSOLUTE REPO R CODE

def calculate_mutation_beta_densities(mutation_df,
                                    alt_count_col='alt', 
                                    ref_count_col='ref',
                                    # alt_count_col='t_alt_count', 
                                    # ref_count_col='t_ref_count',
                                    n_grid=100):
    """
    
    """
    # Initialize the mutation grid
    mut_grid = np.zeros((mutation_df.shape[0], n_grid))
    
    # Compute the coverage and allele frequency
    # MIGHT CHANGE COVERAGE VARIABLE NAME TO => TOTAL READ COUNT
    coverage = mutation_df[alt_count_col] + mutation_df[ref_count_col]
    allele_fraction = mutation_df[alt_count_col] / coverage
    
    # Generate the grid values
    grid_vals = np.arange(1, n_grid + 1) / (n_grid + 1)
    
    # Calculate the Beta densities for each /sample
    for i in range(mutation_df.shape[0]):
        
        mut_grid[i, :] = beta.pdf(grid_vals, 
                                  mutation_df[alt_count_col][i] + 1,
                                  mutation_df[ref_count_col][i] + 1 
                                #   coverage[i] - mutation_df[alt_count_col][i]+ 1 
                                  )/ n_grid
                            # beta.pdf(grid_vals, 
                            #       coverage[i] * allele_fraction[i] + 1, 
                            #       coverage[i] * (1 - allele_fraction[i] + 1), 
                            #       ) / n_grid
    
    # Check for NaN values (rows with NaN values)
    bad_rows = np.any(np.isnan(mut_grid), axis=1)

    if np.any(bad_rows):
        raise ValueError("NaN values detected in the grid")
    
    return mut_grid


def draw_mut_beta_densities(beta_grid, probability_clonal, hz_del_flag, cols, 
                            draw_indv=True, draw_total=True):
    """
    
    """
    n_grid = beta_grid.shape[1]
    grid_vals = np.arange(1, n_grid + 1) / (n_grid + 1)
    probability_subclonal = 1 - probability_clonal
    probability_subclonal[probability_subclonal < 0] = 0  # Handling round-off error
    
    ## Create color palette using Plotly Express
    # color_palette = px.colors.sample_colorscale('RdBu', np.linspace(0, 1, 1000))  # Custom color scale
    # col_scale = len(color_palette)
    # color_range = (0, 1)
    # color_by = probability_clonal

    # Map probability_clonal values to color palette indices
    # palette_idx = np.floor((color_by - color_range[0]) / (color_range[1] - color_range[0]) * (col_scale - 1)).astype(int)
    # mut_colors = np.array(color_palette)[palette_idx]
    
    # Set color to navy for hz_del_flag
    # for i in range(len(hz_del_flag)):
    #     if hz_del_flag[i]:
    #         mut_colors[i] = 'navy'

    # Initialize figure
    fig = go.Figure()
    
    # # Plot individual densities if draw_indv is True
    # if draw_indv:
    #     for i in range(beta_grid.shape[0]):
    #         fig.add_trace(go.Scatter(
    #             x=grid_vals, 
    #             y=beta_grid[i, :], 
    #             mode='lines', 
    #             # CHANGE SO THAT IT IS JUST BLEU COLOR
    #             # line=dict(color=mut_colors[i]), 
    #             line=dict(color="blue"), 
    #             name=f'Sample {i+1}'  # Optional: Add label for each line
    #         ))
    
    # initializes clonal and subclonal weighted grids
    clonal_grid = np.zeros_like(beta_grid)
    subclonal_grid = np.zeros_like(beta_grid)
    # weighted_sum_clonal_grid = np.zeros(len(clonal_grid))
    # weighted_sum_subclonal_grid = np.zeros(len(subclonal_grid))

    for i in range(beta_grid.shape[0]):
        clonal_grid[i, :] = beta_grid[i, :] * probability_clonal[i]
        subclonal_grid[i, :] = beta_grid[i, :] * probability_subclonal[i] #(1 - probability_clonal[i])

        # weighted_sum_clonal_grid[i] = np.sum(clonal_grid[i, :]) / len(clonal_grid[i, :])
        # weighted_sum_subclonal_grid[i] = np.sum(subclonal_grid[i, :]) / len(subclonal_grid[i, :])

    # Plot individual densities if draw_indv is True
    # if draw_indv:
         
    for i in range(beta_grid.shape[0]):
        fig.add_trace(go.Scatter(
            x=np.linspace(0, 1, len(clonal_grid)), 
            y=clonal_grid[i, :], 
            mode='lines', 
            line=dict(color="blue"), 
            name=f'Clonal Sample {i+1}'  # Optional: Add label for each line
        ))

        fig.add_trace(go.Scatter(
            x=np.linspace(0, 1, len(subclonal_grid)), 
            y=subclonal_grid[i, :], 
            mode='lines', 
            # CHANGE SO THAT IT IS JUST BLEU COLOR
            # line=dict(color=mut_colors[i]), 
            line=dict(color="grey"), 
            name=f'Subclonal Sample {i+1}'  # Optional: Add label for each line
        ))
    # debugging_value = clonal_grid
    debugging_value = np.sum(clonal_grid, axis=0) / np.max(np.sum(clonal_grid, axis=0))
    # debugging_value = beta_grid
    # debugging_value = grid_vals
    
    # Plot the total weighted densities if draw_total is True
    if True:
        fig.add_trace(go.Scatter(
            x= np.linspace(0, 1, subclonal_grid.shape[0]),#grid_vals, 
            # x= grid_vals,
            # x= beta_grid, 
            y= np.sum(subclonal_grid, axis=0) / np.max(np.sum(subclonal_grid, axis=0)),
            # y= np.sum(subclonal_grid, axis=1) / np.max(np.sum(subclonal_grid, axis=1)),
            mode='lines', 
            line=dict(color="purple", dash='dash'),  # Dashed line for subclonal
            name='Subclonal Weighted'
        ))
        fig.add_trace(go.Scatter(
            x= np.linspace(0, 1, clonal_grid.shape[0]),#grid_vals, 
            # x= grid_vals, 
            y=np.sum(clonal_grid, axis=0) / np.max(np.sum(clonal_grid, axis=0)),
            # y=np.sum(clonal_grid, axis=1) / np.max(np.sum(clonal_grid, axis=1)),
            mode='lines', 
            line=dict(color="black", dash='dash'),  # Dashed line for clonal
            name='Clonal Weighted'
        ))
        
        # fig.add_trace(go.Scatter(
        #     x= np.linspace(0, 1, len(clonal_grid)),#grid_vals, 
        #     y=np.ones(len(clonal_grid)),
        #     mode='lines', 
        #     line=dict(color="black", dash='dash'),  # Dashed line for clonal
        #     name='Clonal Weighted'
        # ))

        # fig.add_trace(go.Scatter(
        #     x= np.linspace(0, 1, len(subclonal_grid)),#grid_vals, 
        #     y= np.ones(len(subclonal_grid)) * 0.75,
        #     mode='lines', 
        #     line=dict(color="purple", dash='dash'),  # Dashed line for subclonal
        #     name='Subclonal Weighted'
        # ))
    
    # Update layout with titles and axis labels
    fig.update_layout(
        title="Beta Densities",
        xaxis_title="Grid Values",
        yaxis_title="Density",
        template="plotly_white"
    )

    # Show the plot
    fig.show()

    return fig, debugging_value

def mut_af_plot(mut_dat, 
                ssnv_colors, 
                draw_indv,
                mode_color="blue", 
                ):
    """
    
    """
    # DEBUGGING!!
    # print("inside the mut_af_plot function")
    # print("these are the columnn names of the mut_dat dataframe")
    # print(list(mut_dat.columns))
    # print()
    # END OF DEBUGGING!!

    coverage = mut_dat["alt"] + mut_dat["ref"]

    probability_clonal = mut_dat["Pr_somatic_clonal"]
    pr_cryptic_scna = mut_dat["Pr_cryptic_SCNA"]
    ssnv_skew = mut_dat.iloc[0]["SSNV_skew"]  # Assuming this is a scalar

    n_grid = 300
    allele_frac_post_probability = calculate_mutation_beta_densities(mut_dat, n_grid=n_grid)

    # DEBUGGING
    # print("allele_frac_post_probability")
    # print(allele_frac_post_probability)
    # END OF DEBUGGING!!

    grid_vals = np.arange(1, n_grid + 1) / (n_grid + 1)
    grid_mat = np.tile(grid_vals, (allele_frac_post_probability.shape[0], 1))

    # Create the plot using plotly
    fig = go.Figure()

    hz_del_flag = mut_dat["q_hat"] == 0
    fig, debugging_value = draw_mut_beta_densities(allele_frac_post_probability, probability_clonal, hz_del_flag, ssnv_colors, 
                                   draw_indv=draw_indv, draw_total=True)
    # return fig, ()

    # Plot vertical lines and annotations for alpha and SSNV skew
    alpha = mut_dat.iloc[0]["purity"]

    fig.add_trace(go.Scatter(
        x=[alpha / 2, alpha / 2],
        y=[0, 1],
        mode="lines",
        line=dict(width=0.5, dash="dash", color=mode_color),
        name=f'${alpha}_hat / 2$'
    ))

    # Annotate alpha/2
    fig.add_annotation(
        x=alpha / 2,
        y=1.05,
        text="$alpha_hat / 2$",
        showarrow=False,
        font=dict(color=mode_color)
    )

    # Plot the line for SSNV skew * alpha / 2
    fig.add_trace(go.Scatter(
        x=[alpha * ssnv_skew / 2, alpha * ssnv_skew / 2],
        y=[0, 1],
        mode="lines",
        line=dict(width=0.5, dash="dash", color="black"),
        name='f_s_hat dot alpha / 2'
    ))

    # Annotate SSNV skew line
    fig.add_annotation(
        x=alpha * ssnv_skew / 2,
        y=1.05,
        text="f_s_hat dot alpha_hat / 2",
        showarrow=False,
        font=dict(color="black")
    )

    # Update layout with titles and axis labels
    fig.update_layout(
        title="Mutation Allele Frequency (MAF) Plot",
        xaxis_title="Fraction of Alternate Reads",
        yaxis_title="Density",
        template="plotly_white",
        showlegend=True
    )

    # Show the plot
    fig.show()

    return fig, {"af_post_pr": allele_frac_post_probability, "grid_mat": grid_mat}, debugging_value

def get_ssnv_on_clonal_cn_multiplicity_densities(
                                                #  seg_dat, 
                                                 mut_dat, 
                                                 af_post_pr, 
                                                 grid_mat, 
                                                #  verbose=False
                                                 ):
    # Remove mutations on HZdels
    hz_del_flag = mut_dat["q_hat"] == 0
    nix = hz_del_flag

    # Calculate somatic delta
    alpha = mut_dat.iloc[0]["purity"]
    Q = mut_dat["q_hat"]
    som_delta = alpha / (2 * (1 - alpha) + alpha * Q)

    # Initialize matrices for multiplicity
    mult_grid = np.full_like(grid_mat, np.nan, dtype=float)
    mult_dens = np.full_like(af_post_pr, np.nan, dtype=float)

    # Update multiplicity grids for non-HZ del mutations
    mult_grid[~nix, :] = (grid_mat[~nix, :] / som_delta)
    mult_dens[~nix, :] = (af_post_pr[~nix, :] * som_delta)

    return {"mult_dens": mult_dens, "mult_grid": mult_grid}

def multiplicity_plot(
                    #   seg_dat, 
                      mut_dat, 
                      af_post_pr, 
                      grid_mat, 
                    #   SSNV_cols, 
                      mode_color, 
                      draw_indv, 
                      verbose=False):
    
    print("I AM IN THE MULTIPLICITY SSNV PLOT!!")
    
    # Filtering data based on conditions (hz.del.ix, SC_CN.ix, nix)
    hz_del_ix = mut_dat["q_hat"] == 0
    SC_CN_ix = mut_dat["H1"].notna()  # SSNVs on subclonal SCNAs
    nix = hz_del_ix | SC_CN_ix | (mut_dat["alt"] == 0)  # Drop force-called mutations
    
    mut_dat = mut_dat[~nix]
    af_post_pr = af_post_pr[~nix, :]
    grid_mat = grid_mat[~nix, :]

    if len(mut_dat) == 0:
        if verbose:
            print("No valid SSNVs to plot")
        return

    # Get SSNV on clonal CN multiplicity densities
    res = get_ssnv_on_clonal_cn_multiplicity_densities(
                                                    #   seg_dat, 
                                                       mut_dat, af_post_pr, grid_mat, 
                                                    #    verbose=verbose
                                                       )

    mult_dens = res["mult_dens"]
    mult_grid = res["mult_grid"]

    pr_clonal = mut_dat["Pr_somatic_clonal"]
    pr_cryptic_SCNA = mut_dat["Pr_cryptic_SCNA"]
    SSNV_skew = mut_dat.iloc[0]["SSNV_skew"]  # Assuming single value for SSNV_skew

    mult_xlim = 2.5

    # Prepare the data for plotting
    # Create a DataFrame for the individual multiplicity densities
    mult_grid_df = pd.DataFrame(mult_grid, columns=['Multiplicity'])
    for i in range(mult_dens.shape[0]):
        mult_grid_df[f'Individual_{i+1}'] = mult_dens[i, :]

    # Plot using Plotly Express for individual densities
    fig = px.line(mult_grid_df, x='Multiplicity', y=mult_grid_df.columns[1:], 
                  labels={'Multiplicity': 'SSNV Multiplicity', 'value': 'Density'},
                  title="SSNV Multiplicity Densities")

    # Add total density curve if requested
    if draw_indv:
        total_density = np.sum(mult_dens, axis=0)  # Summing densities across individuals
        fig.add_trace(go.Scatter(x=mult_grid, y=total_density, mode='lines', 
                                 name='Total Density', line=dict(color='black', width=2)))

    # Add vertical lines for integer multiplicities (v=1, v=2)
    fig.add_vline(x=1, line=dict(dash='dash', width=2, color=mode_color))
    fig.add_vline(x=2, line=dict(dash='dash', width=2, color=mode_color))

    # Add vertical lines for skew multiplicities (SSNV_skew*1, SSNV_skew*2)
    fig.add_vline(x=SSNV_skew * 1, line=dict(dash='dash', width=2, color='black'))
    fig.add_vline(x=SSNV_skew * 2, line=dict(dash='dash', width=2, color='black'))

    # Update layout for better presentation
    fig.update_layout(
        xaxis_range=[0, mult_xlim],
        yaxis_range=[0, 1],  # Assuming density values are normalized to 1
        xaxis_title="SSNV Multiplicity",
        yaxis_title="Density",
        legend_title="Density Curves"
    )

    fig.show()

    return fig

def classify_somatic_variants(prs, pr_thresh):
    # Create boolean arrays based on threshold condition
    subclonal_ix = prs["Pr_subclonal"] > pr_thresh
    subclonal_wt0_ix = prs["Pr_subclonal_wt0"] > pr_thresh
    clonal_ix = prs["Pr_somatic_clonal"] > pr_thresh
    wt0_ix = prs["Pr_wt0"] > pr_thresh
    ge2_ix = prs["Pr_ge2"] > pr_thresh

    # Additional classifications based on conditions
    clonal_het_ix = prs["Pr_somatic_clonal"] * (1 - prs["Pr_wt0"]) > pr_thresh
    homozygous_ix = prs[["Pr_wt0", "Pr_subclonal_wt0"]].sum(axis=1) > pr_thresh

    # Combine results into a DataFrame
    result = pd.DataFrame({
        "subclonal_ix": subclonal_ix,
        "subclonal_wt0_ix": subclonal_wt0_ix,
        "clonal_ix": clonal_ix,
        "wt0_ix": wt0_ix,
        "clonal_het_ix": clonal_het_ix,
        "ge2_ix": ge2_ix,
        "homozygous_ix": homozygous_ix
    })

    return result

# Define the function for drawing mutation multiplicity densities using plotly
def draw_mut_multiplicity_densities(mut_pr, grid, pr_clonal, pr_cryptic_SCNA, x_lim, xlab, cols, draw_indv=True, draw_total=True):
    def get_grid_combined_mut_densities(mut_pr, pr_clonal, grid, x_lim):
        bin_w = x_lim / 100
        breaks = np.arange(0, x_lim + bin_w, bin_w)
        mult_grid = breaks
        
        # pr-weighted
        grid_dens = np.zeros((mut_pr.shape[0], len(mult_grid)))
        
        for i in range(mut_pr.shape[0]):
            x = grid[i, :]
            y = mut_pr[i, :]
            if np.sum(~np.isnan(y)) > 2:
                # Interpolate the data to match the grid
                interp_func = interp1d(x, y, kind='linear', bounds_error=False, fill_value=np.nan)
                grid_dens[i, :] = interp_func(mult_grid)
        
        y_lim = np.nansum(grid_dens, axis=0)
        return grid_dens, mult_grid, y_lim

    # Subclonal calculation
    pr_subclonal = np.clip(1 - pr_clonal, 0, None)
    
    # Define color palette based on input
    colpal = px.colors.sequential.Plasma  # Using Plotly's color palette
    color_range = (0, 1)
    pal_idx = np.floor(((pr_clonal - color_range[0]) / (color_range[1] - color_range[0])) * (len(colpal) - 1)).astype(int)
    mut_colors = [colpal[idx] for idx in pal_idx]
    
    # Handle cryptic SCNA
    ix = pr_cryptic_SCNA > 0.5
    for i in ix:
        mut_colors[i] = 'mediumorchid2'
    
    # Get combined densities for clonal and subclonal cases
    clonal_dens, mult_grid, _ = get_grid_combined_mut_densities(mut_pr, pr_clonal, grid, x_lim)
    sc_dens, _, _ = get_grid_combined_mut_densities(mut_pr, pr_subclonal, grid, x_lim)

    fig = go.Figure()

    # Draw individual mutation densities if required
    if draw_indv:
        for i in range(mut_pr.shape[0]):
            fig.add_trace(go.Scatter(x=mult_grid, y=clonal_dens[i, :], mode='lines', line=dict(color=mut_colors[i])))

    # Set NAs to 0 for clonal and subclonal densities
    clonal_dens = np.nan_to_num(clonal_dens)
    sc_dens = np.nan_to_num(sc_dens)

    # Draw total clonal and subclonal densities
    if draw_total:
        ncl = np.sum(clonal_dens * pr_clonal, axis=0)
        nsbcl = np.sum(sc_dens * pr_subclonal, axis=0)

        fig.add_trace(go.Scatter(x=mult_grid, y=ncl / np.max(ncl), mode='lines', line=dict(color=cols[1], dash='dash'), name='Clonal'))
        fig.add_trace(go.Scatter(x=mult_grid, y=nsbcl / np.max(nsbcl), mode='lines', line=dict(color=cols[0], dash='dash'), name='Subclonal'))

    fig.update_layout(
        title="Mutation Multiplicity Densities",
        xaxis_title=xlab,
        yaxis_title="Density",
        xaxis=dict(range=[0, x_lim]),
        showlegend=True
    )

    return fig