import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
from cnv_suite.visualize import plot_acr_interactive, plot_acr_subplots
import time

csize = {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260,
        '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747,
        '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
        '16': 90354753, '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
        '21': 48129895, '22': 51304566, '23': 156040895, '24': 57227415}

# cmesser: https://github.com/getzlab/cnv_suite/blob/f88d0bc285a880c2805553762ab939f12e662ad6/cnv_suite/utils/cnv_helper_methods.py
def calc_cn_levels(purity, ploidy, avg_cn=1):
    """Calculate CN zero line and difference between CN levels based on given purity, ploidy and average.
    :param purity: sample tumor purity
    :param ploidy: sample tumor ploidy
    :param avg_cn: average CN value across genome, default 1
    :return: (CN_zero_value, CN_delta_value)
    """
    avg_ploidy = purity * ploidy + 2 * (1 - purity)
    cn_delta = avg_cn * 2 * purity / avg_ploidy
    cn_zero = avg_cn * 2 * (1 - purity) / avg_ploidy
    return cn_zero, cn_delta


def plot_cnp_histogram(
    seg_df,
    mu_major_col,
    mu_minor_col,
    length_col,
    max_mu=2,
    step=0.05,
    fig=None,
    fig_row=None,
    fig_col=None,
):
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
      
def gen_mut_figure(maf_fn,
                   chromosome_col='Chromosome', 
                   start_position_col='Start_position', 
                   hugo_symbol_col='Hugo_Symbol',
                   variant_type_col='Variant_Type',
                   alt_count_col='t_alt_count',
                   ref_count_col='t_ref_count',
                   hover_data=[]  # TODO: include
                  ):
    fig = make_subplots(rows=1, cols=1)
    maf_df = pd.read_csv(maf_fn, sep='\t', encoding='iso-8859-1')
    if maf_df[chromosome_col].dtype == 'object':
        maf_df[chromosome_col].replace({'X': 23, 'Y': 24}, inplace=True)
    maf_df[chromosome_col] = maf_df[chromosome_col].astype(str)
    
    maf_df['new_position'] = maf_df.apply(lambda r: csize[r[chromosome_col]] + r[start_position_col], axis=1)
    maf_df['tumor_f'] = maf_df[alt_count_col] / (maf_df[alt_count_col] + maf_df[ref_count_col])
    
    # color by clonal/subclonal
    if len(hover_data) > 0:
        fig = px.scatter(maf_df, x='new_position', y='tumor_f', marginal_y='histogram', hover_data=hover_data)
    else:
        fig = px.scatter(maf_df, x='new_position', y='tumor_f', marginal_y='histogram')
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
    fig.update_yaxes(range=[0, 1])
    return fig

def gen_cnp_figure(acs_fn,
                   sigmas=True, 
                   mu_major_col='mu.major', 
                   mu_minor_col='mu.minor', 
                   length_col='length',
#                    csize=csize
                  ):
    
    seg_df = pd.read_csv(acs_fn, sep='\t', encoding='iso-8859-1')

    acr_fig, _, _, _ = plot_acr_interactive(seg_df, csize, sigmas=sigmas)
    
    hist_fig = plot_cnp_histogram(
        seg_df,
        mu_major_col,
        mu_minor_col,
        length_col
    )

    cnp_fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.77, 0.25])
    for t in acr_fig.data:
        cnp_fig.add_trace(t, row=1, col=1)

    for t in hist_fig.data:
        cnp_fig.add_trace(t, row=1, col=2)

    cnp_fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")

    return cnp_fig


def parse_absolute_soln(rdata_path: str): # has to be a local path   
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
    end = time.time()
    return mod_tab_df
    
def validate_purity(x):
    return (x >=0) and (x <= 1)

def validate_ploidy(x):
    return (x >=0)
