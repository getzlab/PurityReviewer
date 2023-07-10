import pandas as pd
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import plotly.graph_objects as go

from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from cnv_suite.visualize import plot_acr_interactive
from PurityReviewers.AppComponents.utils import gen_cnp_figure, gen_mut_figure, CSIZE_DEFAULT, parse_absolute_soln


absolute_rdata_cols = ['alpha', 'tau', 'tau_hat', '0_line', '1_line',
                       'sigma_H', 
                       'theta_Q', 
                       'lambda',  
                       'SCNA_likelihood', 
                       'Kar_likelihood', 
                       'SSNVs_likelihood']


def gen_absolute_solutions_report_new_data(
    data: GenericData,
    data_id, 
    selected_row_array, # dash app parameters come first
    rdata_fn_col,
    acs_col, 
    maf_col,
    mut_fig_hover_data,
    csize=None,
    custom_parse_absolute_soln=None
):
    
    data_df = data.df
    r = data_df.loc[data_id]
    
    parse_absolute_soln_func = custom_parse_absolute_soln if custom_parse_absolute_soln is not None else parse_absolute_soln
    try:
        absolute_rdata_df = parse_absolute_soln_func(r[rdata_fn_col])
    except Exception as e:
        print(e)
        absolute_rdata_df = pd.DataFrame()

    absolute_rdata_df = absolute_rdata_df.round(2)
    
    cnp_fig = gen_cnp_figure(r[acs_col], csize=CSIZE_DEFAULT)
    mut_fig = gen_mut_figure(r[maf_col], hover_data=mut_fig_hover_data, csize=CSIZE_DEFAULT)

    # add 1 and 0 lines
    mut_fig_with_lines = go.Figure(mut_fig)
    cnp_fig_with_lines = go.Figure(cnp_fig)
    
    purity = 0
    ploidy = 0
    
    if absolute_rdata_df.shape[0] > 0:
        solution_data = absolute_rdata_df.iloc[selected_row_array[0]]
        i = 0
        line_height = solution_data['0_line']
        while line_height < 2:
            line_height = solution_data['0_line'] + (solution_data['step_size'] * i)
            cnp_fig_with_lines.add_hline(y=line_height, 
                                         line_dash="dash", 
                                         line_color='black',
                                         line_width=1
                                        )
            i += 1

        half_1_line = solution_data['alpha'] / 2.0
        mut_fig_with_lines.add_hline(y=half_1_line, 
                                    line_dash="dash", 
                                    line_color='black',
                                    line_width=1)

        mut_fig_with_lines.update_yaxes(range=[0, half_1_line * 2])
        
        purity = solution_data['alpha']
        ploidy = solution_data['tau_hat']

    return [absolute_rdata_df.to_dict('records'), 
            cnp_fig_with_lines, 
            mut_fig_with_lines,
            purity,
            ploidy, 
            [0],
            0
            ]


def gen_absolute_solutions_report_internal(
    data: GenericData,
    data_id,
    selected_row_array,  # dash app parameters come first
    rdata_fn_col,
    acs_col, 
    maf_col,
    mut_fig_hover_data,
    csize=None,
    custom_parse_absolute_soln=None,
):
    output_data = gen_absolute_solutions_report_new_data(
        data,
        data_id,
        selected_row_array,  # dash app parameters come first
        rdata_fn_col,
        acs_col, 
        maf_col,
        mut_fig_hover_data,
        csize=CSIZE_DEFAULT,
        custom_parse_absolute_soln=custom_parse_absolute_soln,
    )
    output_data[-2] = selected_row_array
    output_data[-1] = selected_row_array[0] + 1 # 1 indexed
    return output_data


def gen_absolute_solutions_report_layout():
    return html.Div(
        children=[
            html.H2('Absolute Solutions Table'),
            dash_table.DataTable(
               id='absolute-rdata-select-table',
               columns=[
                   {"name": i,
                    "id": i} for i in absolute_rdata_cols
               ],
               data=pd.DataFrame(columns=absolute_rdata_cols).to_dict(
                   'records'),
               editable=False,
               filter_action="native",
               sort_action="native",
               sort_mode="multi",
               row_selectable="single",
               row_deletable=False,
               selected_columns=[],
               selected_rows=[0],
               page_action="native",
               page_current=0,
               page_size=5),
            html.H2('Copy Number Profile'),
            html.Div([html.P('Copy number solution: ',
                             style={'display': 'inline'}),
                      html.P(0, id='absolute-solution-idx',
                             style={'display': 'inline'})]),
            html.Div([html.P('Purity: ',
                            style={'display': 'inline'}),
                     html.P(0, id='absolute-purity',
                            style={'display': 'inline'})]),
            html.Div([html.P('Ploidy: ',
                            style={'display': 'inline'}),
                     html.P(0, id='absolute-ploidy',
                            style={'display': 'inline'})]),
            dcc.Graph(id='cnp-graph', figure={}),
            dcc.Graph(id='mut-graph', figure={})
        ]
    )


def gen_absolute_solutions_report_component():
    
    return AppComponent(
        'Absolute Solutions',
        layout=gen_absolute_solutions_report_layout(),
        new_data_callback=gen_absolute_solutions_report_new_data,
        internal_callback=gen_absolute_solutions_report_internal,
        callback_input=[Input('absolute-rdata-select-table', 'selected_rows')],
        callback_output=[
            Output('absolute-rdata-select-table', 'data'),
            Output('cnp-graph', 'figure'),
            Output('mut-graph', 'figure'),
            Output('absolute-purity', 'children'),
            Output('absolute-ploidy', 'children'),
            Output('absolute-rdata-select-table', 'selected_rows'),
            Output('absolute-solution-idx', 'children')
        ],
    )
    
    
    