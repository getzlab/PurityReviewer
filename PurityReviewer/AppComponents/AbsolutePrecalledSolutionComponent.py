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
import pandas as pd

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
    # line_0,
    # line_1,
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

    step_size: float, default=0.05
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
    # defaults to selecting all precalled purity values
    purity_range_lower = 0
    purity_range_upper = 100

    if step_size is None:
        step_size = 0.05

    data_df = data.df
    r = data_df.loc[data_id] # gets a specific row of sample data
    cnp_fig = gen_cnp_figure(r[acs_col], csize=CSIZE_DEFAULT)
    
    if manual_input_source == "Select All":
        
        # displays all precalled purity values!!
        line_0 = 0
        line_1 = 50 

        slider_value = [line_0, line_1]
    else:
        
        if manual_input_source == "Use slider":
            line_0 = slider_value[0]
            line_1 = slider_value[1]

        # gets the lower and upper range of purity values
        purity_range_lower = purity - (line_1 - line_0)
        purity_range_upper = purity + (line_1 - line_0)

        # checks to make sure you dont get negative purity values
        if purity_range_lower < 0:
            purity_range_lower = 0

        # checks to make sure purity value isn't greater than 100%
        if purity_range_upper > 100:
            purity_range_upper = 100
        
    # add 1 and 0 lines
    cnp_fig_with_lines = go.Figure(cnp_fig)
    i = 0
    line_height = line_0
    line_difference = line_1 - line_0

    # creates a horizontal slide
    #MODIFY THAT IT PROPERLY CREATES A HORIZONTAL SLIDER
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
        line_1,
        purity_range_lower,
        purity_range_upper
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
        step_size = 0.05
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
                            ])
                        ]),
                    dbc.Row(
                        [
                          
                            # creating a horizontal slider for selecting a range of purity values
                            dcc.RangeSlider(
                                id='custom-cnp-slider', min=0.0, max=50.0,
                                step=step_size,
                                allowCross=False,
                                value=[0, 5], # initial value on slider 
                                marks={i: f'{i}' for i in range(0, 50, 5)}, 
                                vertical=False, # want to make a horizontal slider
                                tooltip={"placement": "right", "always_visible": True}
                            )
                        ], 
                        md=10
                    )
                ]
            )
        ]

def filter_purity_values(data, purity_val_lower_range, purity_val_upper_range, 
                         current_purity_val, purity_val_col) -> GenericData:
    """
    Only returns the purity values within a specific range

    Parameters:
        data: GenericData, 
        purity_value_lower_range:
        purity_value_upper_range:
        purity_val_col: str, name of the column that has the purity value

    Return:
        returns a generic data object with only rows that contained purity values 
    in a specific range
    """
    data = data.copy()

    data_df = data.df
    annotation_df = data.annot_df
    purity_val_range = purity_val_upper_range - purity_val_lower_range

    filter_row_idx = pd.where((data_df[purity_val_col] >= current_purity_val - purity_val_range) and (data_df[purity_val_col] <= current_purity_val + purity_val_range))

    # only gets the purity values within a specific range of the current purity value
    data.df = data_df.iloc[filter_row_idx]
    data.annot_df = annotation_df.iloc[filter_row_idx]

    return data

def gen_absolute_custom_solution_component(step_size=None):
    """
    Generates an AppComponent defining the interactive elements for setting a manual purity/ploidy solution using a copy number profile

    Returns
    =======
    AnnoMate.AppComponent
        AppComponent defining the interactive elements for setting a manual purity/ploidy solution using a copy number profile
    """
    
    return AppComponent(
        'Precalled Purity',
        layout= gen_precalled_absolute_custom_solution_layout(step_size=step_size),
        new_data_callback=gen_custom_precalled_absolute_component,
        internal_callback=gen_custom_precalled_absolute_component,
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
    
