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

def calculate_multiplicity(maf_df, alpha):
    # Local copy number
    q_values = maf_df['q_hat']

    # calculates Allele fraction
    allele_fraction = maf_df['alt'] / (maf_df['alt']+maf_df['ref'])

    # af = (alpha * multiplicity) / (alpha * q + (1-alpha)*2)
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

def calculate_mutation_beta_densities(
                                    mutation_df,
                                    alt_count_col='alt', 
                                    ref_count_col='ref',
                                    n_grid=100
                                    ):
    """
    Calculate the mutation beta distribution for the mutation read

    Parameters
    ==========
    mutation_df: pd.Dataframe
        mutation dataframe

    alt_count_col: str
        the name of the alternate read count column

    ref_count_col: str
        the name of the reference read count column

    n_grid: int
        the number of columns for the mutation matrix

    Return
    ======
    2D array 
        array of the beta distribution for each sample
    """
    # Initialize the mutation grid
    beta_distribution_matrix = np.zeros((mutation_df.shape[0], n_grid))
    
    total_mutation_reads = mutation_df[alt_count_col] + mutation_df[ref_count_col]    
    
    # Generate the grid values
    normalized_array = np.arange(1, n_grid + 1) / (n_grid + 1) # a normalized array of values from 0 to 1
    
    # Calculate the Beta densities for each sample
    for i in range(mutation_df.shape[0]):
        
        beta_distribution_matrix[i, :] = beta.pdf(normalized_array, 
                                         mutation_df[alt_count_col][i] + 1,
                                         mutation_df[ref_count_col][i] + 1  # num of ref count mutations + 1 = total count of mutations - num of alt count mutations + 1
                                        ) / n_grid # makes the distribution between 0 and 1
    
    # Check for NaN values (rows with NaN values)
    bad_rows = np.any(np.isnan(beta_distribution_matrix), axis=1)

    if np.any(bad_rows):
        raise ValueError("NaN values detected in the grid")
    
    return beta_distribution_matrix

def draw_mut_beta_densities(
                            beta_grid, 
                            probability_clonal, 
                            ):
    """
    Creates the figure with the individual clonal and subclonal beta distribution lines, and the 
        weight clonal and subclonal dashed lines for the fraction alternate read plot and returns it

    Parameters
    ==========
    beta_grid: array
        2D array of the allele fraction posterior probability
    
    probability_clonal: array
        1D array of the probability of a given sample being clonal
    
    Return
    ======
    go.Figure
        a figure with the individual clonal and subclonal lines and the dashed clonal and subclonal lines
    """
    beta_grid = beta_grid.copy()

    indices_to_remove = np.isnan(probability_clonal)
    probability_clonal_cleaned = probability_clonal[~indices_to_remove]
    
    probability_subclonal_cleaned = 1 - probability_clonal_cleaned
    # probability_subclonal = 1 - probability_clonal
    probability_subclonal_cleaned[probability_subclonal_cleaned < 0] = 0  # Handling round-off error
    beta_grid_cleaned = beta_grid[~indices_to_remove] # removes nan values from the beta distribution matrix
    
    # Initialize figure
    fig = go.Figure()
    
    # initializes clonal and subclonal weighted grids
    clonal_grid = np.zeros_like(beta_grid)
    subclonal_grid = np.zeros_like(beta_grid)

    # TRY TO IMPROVE THE SPEED OF THIS BY USING * INSTEAD OF FOR LOOP!!!
    # probability_clonal_vector = probability_clonal_cleaned.to_numpy()
    # probability_subclonal_vector = probability_subclonal_cleaned.to_numpy()
    # probability_subclonal_vector = np.reshape(probability_subclonal_vector, (probability_subclonal_vector.shape[0], 1))
    # probability_clonal_vector = np.reshape(probability_clonal_vector, (probability_clonal_vector.shape[0], 1))
    
    # subclonal_grid = beta_grid_cleaned * probability_subclonal_vector
    # clonal_grid  = beta_grid_cleaned * probability_clonal_vector

    for i in range(beta_grid_cleaned.shape[0]):
        clonal_grid[i, :] = beta_grid_cleaned[i, :] * probability_clonal_cleaned[i]
        subclonal_grid[i, :] = beta_grid_cleaned[i, :] * probability_subclonal_cleaned[i] 

    # Plot individual densities
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
            line=dict(color="grey"), 
            name=f'Subclonal Sample {i+1}'  # Optional: Add label for each line
        ))
    
    # Plots the total weighted densities
    fig.add_trace(go.Scatter(
        x= np.linspace(0, 1, subclonal_grid.shape[0]),
        y= np.sum(subclonal_grid, axis=0) / np.max(np.sum(subclonal_grid, axis=0)),
        mode='lines', 
        line=dict(color="grey", dash='dash'),  # Dashed line for subclonal
        name='Subclonal Weighted'
    ))
    fig.add_trace(go.Scatter(
        x= np.linspace(0, 1, clonal_grid.shape[0]),
        y=np.sum(clonal_grid, axis=0) / np.max(np.sum(clonal_grid, axis=0)),
        mode='lines', 
        line=dict(color="blue", dash='dash'),  # Dashed line for clonal
        name='Clonal Weighted'
    ))
      
    # Show the plot
    fig.show()

    return fig

def mut_allele_fraction_plot(
                             mutation_df, 
                            ):
    """
    Generates the muation allele fraction plot and returns the figure for the plot

    Parameters
    ==========
    mutation_df: pd.Dataframe
        the data

    Return
    ======
        go.Figure
    
    """
    # coverage = mutation_df["alt"] + mutation_df["ref"]

    # gets the probabilities for subclonal and clonal alternate reads
    mode_color="blue"
    probability_clonal = mutation_df["Pr_somatic_clonal"]
    # pr_cryptic_scna = mutation_df["Pr_cryptic_SCNA"]

    # gets ssnv skew and purity value
    # ssnv_skew = mut_dat.iloc[0]["SSNV_skew"]  
    # alpha = mut_dat.iloc[0]["purity"] 
    ssnv_skew = np.mean(mutation_df["SSNV_skew"])  
    alpha = np.mean(mutation_df["purity"]) 

    n_grid = 300
    allele_frac_post_probability = calculate_mutation_beta_densities(mutation_df, n_grid=n_grid)
    normalized_array = np.arange(1, n_grid + 1) / (n_grid + 1) # gets a normalized array of with values between 0 and 1
    grid_mat = np.tile(normalized_array, (allele_frac_post_probability.shape[0], 1))

    # Create the plot using plotly
    fig = go.Figure()
    fig = draw_mut_beta_densities(allele_frac_post_probability, 
                                  probability_clonal, 
                                 )
    
    # plots the vertical line for alpha/2 on figure
    fig.add_trace(go.Scatter(
        x=[alpha / 2, alpha / 2],
        y=[0, 1],
        mode="lines",
        line=dict(width=0.5, dash="dash", color=mode_color),
        name='alpha_hat / 2'
    ))

    # Adds alpha/2 label to figure
    fig.add_annotation(
        x=alpha / 2,
        y=1.05,
        text="alpha_hat / 2",
        showarrow=False,
        font=dict(color=mode_color)
    )

    # plots the SSNV skew * alpha / 2 line on the figure
    fig.add_trace(go.Scatter(
        x=[alpha * ssnv_skew / 2, alpha * ssnv_skew / 2],
        y=[0, 1],
        mode="lines",
        line=dict(width=0.5, dash="dash", color="black"),
        name='f_s_hat dot alpha / 2'
    ))

    # Adds SSNV skew line label to figure
    fig.add_annotation(
        x=alpha * ssnv_skew / 2,
        y=1.15,
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
        showlegend=True,
        xaxis=dict(range=[0, 1]),
        yaxis=dict(range=[0, 1.2]) # increase the y limit so the skew and purity labels do not get cut off
    )

    # Show the plot
    fig.update_layout(showlegend=False)
    fig.show()

    return fig, {"af_post_pr": allele_frac_post_probability, "grid_mat": grid_mat}

def get_ssnv_on_clonal_cn_multiplicity_densities(
                                                #  seg_dat, 
                                                 mutation_df, 
                                                 af_post_pr, 
                                                 grid_mat, 
                                                #  verbose=False
                                                 ):
    # Remove mutations on HZdels
    hz_del_flag = mutation_df["q_hat"] == 0
    nix = hz_del_flag

    # Calculate somatic delta
    alpha = mutation_df.iloc[0]["purity"]
    Q = mutation_df["q_hat"]
    som_delta = alpha / (2 * (1 - alpha) + alpha * Q)
    grid_mat_cleaned = grid_mat[~nix, :]
    af_post_pr_cleaned = af_post_pr[~nix, :]

    # Initialize matrices for multiplicity
    multiplicity_grid = np.full_like(grid_mat_cleaned, np.nan, dtype=float)
    multiplicity_dens = np.full_like(af_post_pr_cleaned, np.nan, dtype=float)

    # Update multiplicity grids for non-HZ del mutations
    # multiplicity_grid[~nix, :] = (grid_mat[~nix, :] / som_delta)
    # multiplicity_dens[~nix, :] = (af_post_pr[~nix, :] * som_delta)
    som_delta = som_delta.to_numpy()
    som_delta_vector = np.reshape(som_delta, (som_delta.shape[0], 1))

    multiplicity_grid = (grid_mat_cleaned / som_delta_vector)
    multiplicity_dens = (af_post_pr_cleaned * som_delta_vector)

    return {"mult_dens": multiplicity_dens, "mult_grid": multiplicity_grid}

def gen_multiplicity_plot(
                      mutation_df, 
                      af_post_pr, 
                      grid_mat, 
                    ):    
    """ 
    
    """
    # Filtering data based on conditions (hz.del.ix, SC_CN.ix, alternate read counts of zero)
    hz_del_ix = mutation_df["q_hat"] == 0
    SC_CN_ix = mutation_df["H1"].notna()  # SSNVs on subclonal SCNAs
    nix = hz_del_ix | SC_CN_ix | (mutation_df["alt"] == 0)  # Drop force-called mutations
    
    mutation_df_cleaned = mutation_df[~nix]
    af_post_pr_cleaned = af_post_pr[~nix, :]
    grid_mat_cleaned = grid_mat[~nix, :]

    # checks if there are no valid SSNVs to ploit
    if len(mutation_df_cleaned) == 0:
        # throws an error
        raise Exception("No valid SSNVs to plot")
         
    # Get SSNV on clonal CN multiplicity densities
    ssnv_multiplicities = get_ssnv_on_clonal_cn_multiplicity_densities(
                                                       mutation_df_cleaned, 
                                                       af_post_pr_cleaned, 
                                                       grid_mat_cleaned, 
                                                       )

    mult_dens = ssnv_multiplicities["mult_dens"]
    mult_grid = ssnv_multiplicities["mult_grid"]

    pr_clonal = mutation_df_cleaned["Pr_somatic_clonal"]
    pr_cryptic_SCNA = mutation_df_cleaned["Pr_cryptic_SCNA"]
    SSNV_skew = mutation_df_cleaned.iloc[0]["SSNV_skew"]  # Assuming single value for SSNV_skew

    mult_xlim = 2.5

    fig = go.Figure()

    # # Prepare the data for plotting
    # # Create a DataFrame for the individual multiplicity densities
    # multiplicity_densities_df = pd.DataFrame(mult_grid.to_numpy(), columns=['Multiplicity'])
    
    # draws the individal subclonal and clonal ssnv multiplicity lines!!
    fig = draw_mut_multiplicity_densities(mult_dens, 
                                    mult_grid, 
                                    pr_clonal, 
                                    pr_cryptic_SCNA, 
                                    mult_xlim, 
                                    )

    # Add total density curve if requested
    total_density = np.sum(mult_dens, axis=0) / np.max(np.sum(mult_dens, axis=0))  # normalized sum of densities across individuals
    debugging_value = total_density

    # Add vertical lines for integer multiplicities (v=1, v=2)
    fig.add_trace(go.Scatter(
        x=[1, 1],
        y=[0, 1],
        mode="lines",
        line=dict(width=2, dash="dash", color="green"),
        name="multiplicity line at 1"
    ))

    fig.add_trace(go.Scatter(
        x=[2, 2],
        y=[0, 1],
        mode="lines",
        line=dict(width=2, dash="dash", color="green"),
        name="multiplicity line at 2"
    ))

    fig.add_trace(go.Scatter(
        x=[SSNV_skew * 1, SSNV_skew * 1],
        y=[0, 1],
        mode="lines",
        line=dict(width=2, dash="dash", color="black"),
        name="ssnv_skew * 1"
    ))

    fig.add_trace(go.Scatter(
        x=[SSNV_skew * 2, SSNV_skew * 2],
        y=[0, 1],
        mode="lines",
        line=dict(width=2, dash="dash", color="black"),
        name="ssnv_skew * 2"
    ))

    # Adds SSNV skew line label to figure
    fig.add_annotation(
        x=SSNV_skew * 1,
        y=1.15,
        text="ssnv_skew * 1",
        showarrow=False,
        font=dict(color="black")
    )

    # Adds SSNV skew line label to figure
    fig.add_annotation(
        x=SSNV_skew * 2,
        y=1.15,
        text="ssnv_skew*2",
        showarrow=False,
        font=dict(color="black")
    )

    # Update layout for better presentation
    fig.update_layout(
        xaxis_range=[0, mult_xlim],
        yaxis_range=[0, 1.2],  # Assuming density values are normalized to 1
        xaxis_title="SSNV Multiplicity",
        template="plotly_white",
        yaxis_title="Density",
        legend_title="Density Curves"
    )

    fig.show()
    fig.update_layout(showlegend=False)

    return fig, debugging_value

def draw_mut_multiplicity_densities(mut_pr, 
                                    grid, 
                                    pr_clonal, 
                                    pr_cryptic_SCNA, 
                                    x_lim, 
                                    ):
    
    def get_grid_combined_mut_densities(
                                        mut_pr, 
                                        pr_clonal, 
                                        grid, 
                                        x_lim
                                        ):
        
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
    pr_subclonal = 1 - pr_clonal
    pr_subclonal_vector = np.clip(1 - pr_clonal, 0, None)
    pr_subclonal_vector = pr_subclonal.to_numpy()
    pr_clonal_vector = pr_clonal.to_numpy()
    
    # Get combined densities for clonal and subclonal cases
    clonal_dens, clonal_mult_grid, _ = get_grid_combined_mut_densities(mut_pr, pr_clonal_vector, grid, x_lim)
    subclonal_dens, subclonal_mult_grid, _ = get_grid_combined_mut_densities(mut_pr, pr_subclonal_vector, grid, x_lim)

    pr_clonal_vector = np.reshape(pr_clonal_vector, (pr_clonal_vector.shape[0], 1))
    pr_subclonal_vector = np.reshape(pr_subclonal_vector, (pr_subclonal_vector.shape[0], 1))

    # Set NAs to 0 for probability of being clonal or subclonal
    pr_clonal_vector = np.nan_to_num(pr_clonal_vector)
    pr_subclonal_vector = np.nan_to_num(pr_subclonal_vector)

    # Set NAs to 0 for clonal and subclonal beta distributions
    clonal_dens = np.nan_to_num(clonal_dens)
    sc_dens = np.nan_to_num(subclonal_dens)

    fig = go.Figure()

    # Draw individual subclonal and clonal mutation densities
    for i in range(mut_pr.shape[0]):
        fig.add_trace(go.Scatter(
            # x=clonal_mult_grid,
            x=np.linspace(0, 2.5, len(clonal_dens)), 
            y=clonal_dens[i, :] * pr_clonal_vector[i], 
            mode='lines', 
            line=dict(color="blue"),
            name=f"Clonal Sample {i+1}")
        )
    
        fig.add_trace(go.Scatter(
            # x=subclonal_mult_grid, 
            x=np.linspace(0, 2.5, len(subclonal_dens)),
            y=subclonal_dens[i, :] * pr_subclonal_vector[i], 
            mode='lines', 
            line=dict(color="grey"),
            name=f"Subclonal Sample {i+1}")
        )

    # Draw weighted clonal and subclonal densities
    ncl = np.sum(clonal_dens * pr_clonal_vector, axis=0)
    nsbcl = np.sum(sc_dens * pr_subclonal_vector, axis=0)

    fig.add_trace(go.Scatter(
        # x=clonal_mult_grid,
        x=np.linspace(0, 2.5, len(ncl)), 
        y=ncl / np.max(ncl), 
        mode='lines', 
        line=dict(color="blue", dash='dash'), 
        name='Weighted Clonal')
    )
    fig.add_trace(go.Scatter(
        # x=clonal_mult_grid, 
        x = np.linspace(0, 2.5, len(nsbcl)),
        y=nsbcl / np.max(nsbcl), 
        mode='lines', 
        line=dict(color="grey", dash='dash'), 
        name='Weighted Subclonal')
    )

    fig.update_layout(
        title="Mutation Multiplicity Densities",
        xaxis_title="SSNV multiplicity",
        yaxis_title="Density",
        xaxis=dict(range=[0, x_lim]),
        # showlegend=True
    )

    return fig