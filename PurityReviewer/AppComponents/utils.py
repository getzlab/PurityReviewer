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
from sklearn.preprocessing import MinMaxScaler

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

def Mut_AF_plot(mut_dat, SSNV_cols, mode_color, draw_indv=True):
    # Calculate coverage
    cov = mut_dat['alt'] + mut_dat['ref']

    # Extract the necessary columns
    pr_clonal = mut_dat['Pr_somatic_clonal'].values
    pr_cryptic_SCNA = mut_dat['Pr_cryptic_SCNA'].values
    SSNV_skew = mut_dat['SSNV_skew'].values[0]

    # Generate mutation beta densities (placeholder function)
    n_grid = 300
    af_post_pr = get_mut_beta_densities(mut_dat, n_grid)

    # Create the grid for plotting
    grid = np.linspace(1, n_grid, n_grid) / (n_grid + 1)
    grid_mat = np.tile(grid, (af_post_pr.shape[0], 1))  # Repeat grid for each mutation
    
    # Plot setup
    fig = go.Figure()

    # Draw the densities (individual and/or total)
    hz_del_flag = mut_dat['q_hat'] == 0
    fig = draw_mut_beta_densities(af_post_pr, pr_clonal, hz_del_flag, SSNV_cols, draw_total=True, draw_indv=draw_indv)

    # for trace in traces:
    #     fig.add_trace(trace)

    # Add vertical lines for alpha / 2 and alpha * SSNV_skew / 2
    alpha = mut_dat['purity'].values[0]
    fig.add_shape(
        type='line',
        x0=alpha / 2,
        x1=alpha / 2,
        y0=0,
        y1=1,
        line=dict(color=mode_color, dash='dot', width=2)
    )
    fig.add_annotation(
        x=alpha / 2,
        y=1.05,
        text=r'$\hat{\alpha} / 2$',
        showarrow=False,
        font=dict(color=mode_color)
    )

    fig.add_shape(
        type='line',
        x0=alpha * SSNV_skew / 2,
        x1=alpha * SSNV_skew / 2,
        y0=0,
        y1=1,
        line=dict(color='black', dash='dot', width=2)
    )
    fig.add_annotation(
        x=alpha * SSNV_skew / 2,
        y=1.05,
        text=r'$\hat{f_s} \times \hat{\alpha} / 2$',
        showarrow=False,
        font=dict(color='black')
    )

    # Update layout
    fig.update_layout(
        title="Mutation Allele Fraction Plot",
        yaxis=dict(title="Density"),
        template="plotly_white",
        xaxis_range=[0, 1],
        yaxis_range=[0, 1],
        xaxis=dict(
            title="Fraction of Alternate Reads",
            tickvals=[i/10 for i in range(0, 11, 2)],  # Set the tick positions
            ticktext=[f"{i/10}" for i in range(0, 11, 2)]),  # Set the labels for those tick positions)
        showlegend=False,
        width=400,
    )

    return fig, {"af_post_pr": af_post_pr, "grid_mat": grid_mat}

# Define the function that generates the plot
def draw_mut_beta_densities(beta_grid, pr_clonal, hz_del_flag, cols, draw_indv=True, draw_total=True):
    n_grid = beta_grid.shape[1]
    grid_vals = np.arange(1, n_grid + 1) / (n_grid + 1)  # grid values
    
    pr_subclonal = 1 - pr_clonal
    pr_subclonal[pr_subclonal < 0] = 0  # Fix round-off error
    
    # Create color map based on pr_clonal values
    color_by = pr_clonal
    scaler = MinMaxScaler((0, 1))
    color_by_scaled = scaler.fit_transform(color_by.reshape(-1, 1)).flatten()  # Normalize
    
    # Generate a color palette
    mut_colors = [cols[1] if hz_del else cols[0] for hz_del in hz_del_flag]
    
    # Set up the plot
    fig = go.Figure()

    if draw_indv:
        for i in range(beta_grid.shape[0]):
            fig.add_trace(go.Scatter(
                x=grid_vals,
                y=beta_grid[i, :],
                mode='lines',
                line=dict(color=mut_colors[i]),
                name=f"Individual {i+1}"  # Label each individual line
            ))

    # Pr-weighted clonal and subclonal
    clonal_grid = beta_grid * pr_clonal[:, np.newaxis]
    sc_grid = beta_grid * pr_subclonal[:, np.newaxis]
    
    if draw_total:
        fig.add_trace(go.Scatter(
            x=grid_vals,
            y=np.sum(clonal_grid, axis=0) / np.max(np.sum(clonal_grid, axis=0)),
            mode='lines',
            line=dict(color=cols[1], dash='dash'),
            name="Clonal Total"
        ))
        
        fig.add_trace(go.Scatter(
            x=grid_vals,
            y=np.sum(sc_grid, axis=0) / np.max(np.sum(sc_grid, axis=0)),
            mode='lines',
            line=dict(color=cols[0], dash='dash'),
            name="Subclonal Total"
        ))
    
    return fig

def get_mut_beta_densities(mut_dat, n_grid=100):
    # Initialize mutation grid with zeros
    mut_grid = np.zeros((mut_dat.shape[0], n_grid))
    
    # Calculate coverage and allele frequency
    cov = mut_dat["alt"] + mut_dat["ref"]
    af = mut_dat["alt"] / cov
    
    # Grid values between 0 and 1
    grid_vals = np.linspace(0, 1, n_grid)
    
    # Calculate beta densities for each mutation
    for i in range(mut_dat.shape[0]):
        a = cov[i] * af[i] + 1  # Alpha parameter
        b = cov[i] * (1 - af[i]) + 1  # Beta parameter
        mut_grid[i, :] = beta.pdf(grid_vals, a, b) * 1 / n_grid  # Normalize by n_grid
    
    # Check for NaN values
    if np.any(np.isnan(mut_grid)):
        raise ValueError("NaN values detected in the mutation grid")
    
    return mut_grid

def get_SSNV_on_clonal_CN_multiplicity_densities(seg_dat, mut_dat, af_post_pr, grid_mat, verbose=False):
    # Remove mutations on HZdels
    hz_del_flag = mut_dat["q_hat"] == 0
    nix = hz_del_flag

    alpha = mut_dat.iloc[0]["purity"]  # Assuming purity is stored in the first row
    Q = mut_dat["q_hat"].values
    som_delta = alpha / (2 * (1 - alpha) + alpha * Q)
    # som_delta = som_delta * SSNV_skew^-1  # If SSNV_skew is available, apply this adjustment

    # Initialize empty matrices for mutation grid and densities
    mult_grid = np.full_like(grid_mat, np.nan, dtype=float)
    mult_dens = np.full_like(af_post_pr, np.nan, dtype=float)

    cleaned_grid = grid_mat[~nix, :]
    cleaned_som = som_delta[~nix]
    # cleaned_som_vector = cleaned_som.to_numpy()
    cleaned_som = np.reshape(cleaned_som, (cleaned_som.shape[0], 1))
    cleaned_af_post = af_post_pr[~nix, :]

    # Compute multiplicity grid and densities for non-removed mutations
    mult_grid[~nix, :] = (cleaned_grid / cleaned_som)  # Adjust based on som_delta
    mult_dens[~nix, :] = (cleaned_af_post * cleaned_som)  # Adjust based on som_delta

    return {"mult_dens": mult_dens, "mult_grid": mult_grid}

# Function to draw mutation multiplicity densities
def draw_mut_multiplicity_densities(mult_dens, mult_grid, pr_clonal, pr_cryptic_SCNA, x_lim, cols, xlab, draw_indv, y_lim=1):
    fig = go.Figure()

    # Add individual mutation densities if requested
    if draw_indv:
        for i in range(mult_dens.shape[0]):
            fig.add_trace(go.Scatter(x=mult_grid, y=mult_dens[i, :], mode='lines', line=dict(color=cols[i])))

    # Plot the total density
    total_density = np.nansum(mult_dens * pr_clonal[:, None], axis=0)
    fig.add_trace(go.Scatter(x=mult_grid, y=total_density / np.max(total_density), mode='lines', 
                             line=dict(color=cols[1], dash='dash')))

    # Adjust layout settings
    fig.update_layout(
        title="SSNV Multiplicity Densities",
        # xaxis_title=xlab,
        yaxis_title="Density",
        xaxis=dict(range=[0, x_lim],
                   title="SSNV Multiplicity",
                   tickvals=[i/10 for i in range(0, x_lim*10+1, 5) ],  # Set the tick positions
                   ticktext=[f"{i/10}" for i in range(0, x_lim*10+1, 5)]),  # Set the labels for those tick positions),
        yaxis=dict(range=[0, y_lim]),
        template="plotly_white",
        showlegend=False,
        width=400,
    )

    # Highlight integer multiplicities (1 and 2)
    for v in [1, 2]:
        fig.add_vline(x=v, line=dict(dash="dash", color="black", width=0.5))

    # Highlight skew multiplicities
    SSNV_skew = 1  # Placeholder, use actual skew from data
    for v in [SSNV_skew, SSNV_skew * 2]:
        fig.add_vline(x=v, line=dict(dash="dash", color="red", width=0.5))

    return fig

# Main function to plot multiplicity
def multiplicity_plot(seg_dat, mut_dat, af_post_pr, grid_mat, SSNV_cols, mode_color, draw_indv, verbose=False):
    # Filter out invalid mutations
    hz_del_ix = mut_dat["q_hat"] == 0
    SC_CN_ix = mut_dat["H1"].notna()  # SSNVs on subclonal SCNAs
    nix = hz_del_ix | SC_CN_ix | (mut_dat["alt"] == 0)  # Drop forced-called mutations

    # Apply filtering
    mut_dat_filtered = mut_dat[~nix]
    af_post_pr_filtered = af_post_pr[~nix]
    grid_mat_filtered = grid_mat[~nix]

    if mut_dat_filtered.shape[0] == 0:
        if verbose:
            print("No valid SSNVs to plot")
        return None

    # Get the mutation density data
    res = get_SSNV_on_clonal_CN_multiplicity_densities(seg_dat, mut_dat_filtered, af_post_pr_filtered, grid_mat_filtered, verbose)
    mult_dens = res["mult_dens"]
    mult_grid = res["mult_grid"]

    # Extract required information from the mutation data
    pr_clonal = mut_dat_filtered["Pr_somatic_clonal"].values
    pr_cryptic_SCNA = mut_dat_filtered["Pr_cryptic_SCNA"].values
    SSNV_skew = mut_dat_filtered.iloc[0]["SSNV_skew"]  # Assuming skew is a scalar

    # Plot the densities
    mult_xlim = 2.5
    fig = draw_mut_multiplicity_densities(
        mult_dens, mult_grid, pr_clonal, pr_cryptic_SCNA, x_lim=mult_xlim,
        cols=SSNV_cols, xlab="SSNV multiplicity", draw_indv=draw_indv, y_lim=1
    )
    debugging_value = ""

    return fig, debugging_value

def get_grid_combined_mut_densities(mut_pr, pr_clonal, grid, x_lim):
    bin_w = x_lim / 100
    breaks = np.arange(0, x_lim, bin_w)
    mult_grid = breaks

    grid_dens = np.zeros((mut_pr.shape[0], len(mult_grid)))

    for i in range(mut_pr.shape[0]):
        x = grid[i, :]
        y = mut_pr[i, :]
        if np.sum(~np.isnan(y)) > 2:
            interp_func = interp1d(x, y, kind='linear', bounds_error=False, fill_value=0)
            grid_dens[i, :] = interp_func(mult_grid)

    y_lim = np.nansum(grid_dens, axis=0)
    return grid_dens, mult_grid, y_lim

# Function to draw mutation multiplicity densities using Plotly Dash
def draw_mut_multiplicity_densities(mut_pr, grid, pr_clonal, pr_cryptic_SCNA, x_lim, xlab, cols, 
                                     draw_indv=True, draw_total=True, add=False, y_lim=None):
    # Pr subclonal calculation
    pr_subclonal = 1 - pr_clonal
    pr_subclonal[pr_subclonal < 0] = 0  # Round-off error handling

    # Generate color palette
    color_palette = px.colors.sequential.Blues  # Use the palette from Plotly
    colpal = color_palette[:1000]  # Use the first 1000 colors
    col_scale = len(colpal)
    color_range = [0, 1]
    color_by = pr_clonal
    pal_idx = np.floor((color_by - color_range[0]) / (color_range[1] - color_range[0]) * (col_scale - 1)).astype(int)
    # mut_colors = [colpal[i] for i in pal_idx]

    # Apply color adjustment for cryptic SCNA
    # mut_colors[pr_cryptic_SCNA > 0.5] = "mediumorchid2"

    # Get combined mutation densities for clonal and subclonal
    clonal_dens, mult_grid, _ = get_grid_combined_mut_densities(mut_pr, pr_clonal, grid, x_lim)
    sc_dens, _, _ = get_grid_combined_mut_densities(mut_pr, pr_subclonal, grid, x_lim)

    # Initialize the figure
    fig = go.Figure()

    # If not adding to an existing plot, clear the plot area
    if not add:
        fig.update_layout(
            title="Mutation Multiplicity Densities",
            xaxis_title=xlab,
            yaxis_title="Density",
            xaxis=dict(
                    range=[0, x_lim],
                    tickvals=[i/10 for i in range(0, int(x_lim*10)+4, 5) ],  # Set the tick positions
                    ticktext=[f"{i/10}" for i in range(0,int(x_lim*10)+4, 5)],  # Set the labels for those tick positions),
                    ),
            yaxis=dict(range=[0, y_lim] if y_lim is not None else [0, np.nanmax(clonal_dens)]),
            showlegend=False,
            template="plotly_white",
            width=400,
        )

    # Add individual lines for each mutation
    if draw_indv:
        for i in range(mut_pr.shape[0]):
            fig.add_trace(go.Scatter(x=mult_grid, y=clonal_dens[i, :], mode='lines', line=dict(color="blue")))
                                    #  line=dict(color=mut_colors[i])))

    # Set the densities to 0 where they are NaN
    clonal_dens[np.isnan(clonal_dens)] = 0
    sc_dens[np.isnan(sc_dens)] = 0

    # Add the total densities if requested
    if draw_total:
        cleaned_pr_clonal = np.reshape(pr_clonal, (pr_clonal.shape[0], 1))
        cleaned_pr_subclonal = np.reshape(pr_subclonal, (pr_subclonal.shape[0], 1))

        ncl = np.nansum(clonal_dens * cleaned_pr_clonal, axis=0)
        nsbcl = np.nansum(sc_dens * cleaned_pr_subclonal, axis=0)

        fig.add_trace(go.Scatter(x=mult_grid, y=ncl / np.max(ncl), mode='lines', 
                                 line=dict(color=cols[1], dash='dash')))
        fig.add_trace(go.Scatter(x=mult_grid, y=nsbcl / np.max(nsbcl), mode='lines', 
                                 line=dict(color=cols[0], dash='dash')))

    return fig