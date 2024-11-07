"""
Displays a allelic copy ratio profile with the option to set the 0 and 1 line (via slider or input value) corresponding to an integer assignment to copy number peaks. Automatically calculates the purity given the corresponding solution.
"""
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from cnv_suite.visualize import plot_acr_interactive
from PurityReviewer.AppComponents.utils import gen_cnp_figure, CSIZE_DEFAULT

from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
import os
import pickle
from typing import Union, List, Dict
import sys
from cnv_suite import calc_cn_levels

MANUAL_INPUT_SOURCE = ["Use slider", "Select All"]

# Add a custom component: below, I add a component that allows you to manually set the range of purity values
# NOTE: the purity calculated is tau, NOT tau_g. 
# Use the called purity and tau as inputs to absolute_segforcecall to get tau_g (tau_hat)
def gen_custom_precalled_absolute_component(
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
    """
    Callback function that updates copy number plots when a new sample/pair is selected 
    
    Parameters
    ==========
    data: GenericData
        Data object storing the relevant data to review
        
    data_id: str
        Name of the object being reviewed
        
    slider_value: List[float]
        List of length 2 where the 0-index entry is the current value of the 0-line, and the 1-index entry is the current value of the 1-line

    purity: float
        Current value for purity

    ploidy: float
        Current value for ploidy

    line_0: float
        Current value of the 0-line

    line_1: float
        Current value of the 1-line

    manual_input_source: str ["Manual Purity/ploidy", "Use slider", "Manual 0/1 line"]
        Which mode to set get the purity solution and replot the copy number profile with the corresponding comb
        
    acs_col: str
        Column name in data with path to seg file from alleliccapseg or other tsv with allelic copy ratio measurements.

    step_size: float, default=0.01
        Degree of precision for the slider bar
    
    csize: dict
            Dictionary with chromosome sizes (see AppComponents.utils.CSIZE_DEFAULT for hg19)

    Returns
    =======
    List[float]
        List of length 2 where the 0-index entry is the current/recalculated value of the 0-line, and the 1-index entry is the current value of the 1-line. Used to update the slider bar

    plotly.Figure
        Copy number profile figure with the comb lines plotted on top

    float
        Current/recalculated purity value

    float
        Current/recalculated ploidy value
        
    float
        Current/recalculated 0-line

    float Current/recalculated 1-line
    """
    
    if step_size is None:
        step_size = 0.01
    data_df = data.df
    r = data_df.loc[data_id]
    cnp_fig = gen_cnp_figure(r[acs_col], csize=CSIZE_DEFAULT)
    
    if manual_input_source == "Manual Purity/ploidy":
        cn_0, cn_delta = calc_cn_levels(purity, ploidy)
        
        line_0 = round(cn_0, 2)
        line_1 = round(cn_0 + cn_delta, 2)
        slider_value = [line_0, line_1]
    else:
        
        if manual_input_source == "Use slider":
            line_0 = slider_value[0]
            line_1 = slider_value[1]
            
        # elif manual_input_source == "Manual 0/1 line":
        #     slider_value = [line_0, line_1]
        
        # purity = round((1 - (float(line_0) / float(line_1))) / step_size) * step_size
        # ploidy = round(((2 * (1 - line_0) * (1 - purity)) / (purity * line_0)) / step_size) * step_size  # tau, not tau_g
        
    # add 1 and 0 lines
    cnp_fig_with_lines = go.Figure(cnp_fig)
    i = 0
    line_height = line_0
    line_difference = line_1 - line_0

    # creates a horizontal slide
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
    
def gen_precalled_absolute_custom_solution_layout(step_size=None):
    """
    Generates the layout of the custom purity solutions component in the dashboard

    Returns
    =======
    dash.html
        a plotly dash layout with a copy number plot and multiple options to set the purity and ploidy or set the 0 and 1 line
    """
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
                            #    html.Div([
                            #        html.P('0-line: ', style={'display': 'inline'}),
                            #        html.Div(
                            #            dbc.Input(
                            #                id='custom-cnp-graph-0-line',
                            #                type='number',
                            #                min=0, step=step_size
                            #            ),
                            #            style={'display': 'inline'}
                            #        )
                            #    ]),
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
                                            id='custom-cnp-slider', min=0.0, max=100.0,
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
    """
    Generates an AppComponent defining the interactive elements for setting a manual purity/ploidy solution using a copy number profile

    Returns
    =======
    AnnoMate.AppComponent
        AppComponent defining the interactive elements for setting a manual purity/ploidy solution using a copy number profile
    """
    
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
    