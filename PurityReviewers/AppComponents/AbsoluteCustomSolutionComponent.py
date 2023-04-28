from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from PurityReviewers.AppComponents.utils import gen_cnp_figure
from JupyterReviewer.ReviewDataApp import AppComponent
from JupyterReviewer.DataTypes.GenericData import GenericData
from cnv_suite import calc_cn_levels


MANUAL_INPUT_SOURCE = ["Use slider", "Manual Purity/ploidy", "Manual 0/1 line"]

# Add a custom component: below, I add a component that allows you to manually set the 0 and 1 line combs (cgaitools)
# NOTE: the purity calculated is tau, NOT tau_g. 
# Use the called purity and tau as inputs to absolute_segforcecall to get tau_g (tau_hat)
def gen_custom_absolute_component(
    data: GenericData,
    data_id,
    slider_value,  # dash app parameters come first
    purity,
    ploidy,
    line_0,
    line_1,
    manual_input_source,
    acs_col,
    step_size=None,
    csize=None

):
    if step_size is None:
        step_size = 0.01
    data_df = data.df
    r = data_df.loc[data_id]
    cnp_fig = gen_cnp_figure(r[acs_col], csize=csize)
    
    if manual_input_source == "Manual Purity/ploidy":
        cn_0, cn_delta = calc_cn_levels(purity, ploidy)
        
        line_0 = round(cn_0, 2)
        line_1 = round(cn_0 + cn_delta, 2)
        slider_value = [line_0, line_1]
    else:
        
        if manual_input_source == "Use slider":
            line_0 = slider_value[0]
            line_1 = slider_value[1]
            
        elif manual_input_source == "Manual 0/1 line":
            slider_value = [line_0, line_1]
        
        purity = round((1 - (float(line_0) / float(line_1))) / step_size) * step_size
        ploidy = round(((2 * (1 - line_0) * (1 - purity)) / (purity * line_0)) / step_size) * step_size  # tau, not tau_g
        
    # add 1 and 0 lines
    cnp_fig_with_lines = go.Figure(cnp_fig)
    i = 0
    line_height = line_0
    line_difference = line_1 - line_0
    while line_height < 2:
        line_height = line_0 + (line_difference * i)
        cnp_fig_with_lines.add_hline(y=line_height,
                                     line_dash="dash",
                                     line_color='black',
                                     line_width=1
                                     )
        i += 1

    
    return [
        slider_value,
        cnp_fig_with_lines,
        purity,
        ploidy,
        line_0, 
        line_1
    ]
    
def gen_absolute_custom_solution_layout(step_size=None):
    if step_size is None:
        step_size = 0.01
    return [
            html.Div(
                [
                    dbc.Row([
                        html.Div(
                            [
                                dbc.Label("Choose one"),
                                dbc.RadioItems(
                                    options=[
                                        {
                                            "label": v, 
                                            "value": v
                                        } for v in MANUAL_INPUT_SOURCE
                                    ],
                                    value="Use slider",
                                    id="manual-input-source-radioitems",
                                ),
                            ]
                        )
                    ]),
                    dbc.Row([
                        dbc.Col(
                            [
                                html.Div([
                                    html.P('Purity: ', style={'display': 'inline'}),
                                    html.Div(
                                        dbc.Input(
                                            id='custom-cnp-graph-purity',
                                            type='number',
                                            min=0, max=1.0, step=step_size
                                        ),
                                        style={'display': 'inline'}
                                    )
                                ]),
                                html.Div([
                                    html.P('Ploidy: ', style={'display': 'inline'}),
                                    html.Div(
                                        dbc.Input(
                                            id='custom-cnp-graph-ploidy',
                                            type='number',
                                            min=0, step=step_size
                                        ),
                                        style={'display': 'inline'}
                                    )
                                ]),
                            ],
                           md=2
                        ),
                        dbc.Col(
                            [
                                html.Div([
                                   html.P('1-line: ', style={'display': 'inline'}),
                                   html.Div(
                                       dbc.Input(
                                           id='custom-cnp-graph-1-line',
                                           type='number',
                                           min=0, step=step_size
                                       ),
                                       style={'display': 'inline'}
                                   )
                               ]),
                               html.Div([
                                   html.P('0-line: ', style={'display': 'inline'}),
                                   html.Div(
                                       dbc.Input(
                                           id='custom-cnp-graph-0-line',
                                           type='number',
                                           min=0, step=step_size
                                       ),
                                       style={'display': 'inline'}
                                   )
                               ]),
                            ],
                            md=2
                        )
                    ]),
                    dbc.Row(
                        [
                            dbc.Col(
                                [dcc.Graph(id='custom-cnp-graph', figure={})],
                                md=10
                            ),
                            dbc.Col(
                                [
                                    dbc.Row(
                                        dcc.RangeSlider(
                                            id='custom-cnp-slider', min=0.0, max=2.0,
                                            step=step_size,
                                            allowCross=False,
                                            value=[0.5, 1.0],
                                            marks={i: f'{i}' for i in range(0, 3, 1)},
                                            vertical=True,
                                            tooltip={"placement": "right", "always_visible": True}
                                        )
                                    )
                                ], 
                                md=1
                            )
                        ]
                    ),
                ]
            )
        ]

def gen_absolute_custom_solution_component(step_size=None):
    # Adding another component to prebuilt dash board
    
    return AppComponent(
        'Manual Purity',
        layout=gen_absolute_custom_solution_layout(step_size=step_size),
        new_data_callback=gen_custom_absolute_component,
        internal_callback=gen_custom_absolute_component,
        callback_input=[
            Input('custom-cnp-slider', 'value'),
            Input('custom-cnp-graph-purity', 'value'),
            Input('custom-cnp-graph-ploidy', 'value'),
            Input('custom-cnp-graph-0-line', 'value'),
            Input('custom-cnp-graph-1-line', 'value'),
            Input('manual-input-source-radioitems', 'value')
        ],
        callback_output=[
            Output('custom-cnp-slider', 'value'),
            Output('custom-cnp-graph', 'figure'),
            Output('custom-cnp-graph-purity', 'value'),
            Output('custom-cnp-graph-ploidy', 'value'),
            Output('custom-cnp-graph-0-line', 'value'),
            Output('custom-cnp-graph-1-line', 'value')
        ],
    )
    