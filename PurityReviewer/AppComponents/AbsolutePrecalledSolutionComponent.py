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
    # purity, # Would need to add Input object to callback_inputs in the gen_precalled_custom_component function
    manual_input_source,
    purity_col='PCA__ABSOLUTE__Cancer_DNA_fraction',
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

    line_0: float
        Current value of the 0-line

    line_1: float
        Current value of the 1-line

    manual_input_source: str ["Use slider", "Select All"]
        Which mode to set get the purity solution and replot the copy number profile with the corresponding comb
        
    purity_col: str
        Column name in data with path to precalled purity values

    step_size: float, default=0.05
        Degree of precision for the slider bar
    
    csize: dict
            Dictionary with chromosome sizes (see AppComponents.utils.CSIZE_DEFAULT for hg19)

    Returns
    =======
    List[float]
        List of length 2 where the 0-index entry is the current/recalculated value of the 0-line, and the 1-index entry is the current value of the 1-line. Used to update the slider bar

    float
        Current/recalculated purity value
        
    float
        Current/recalculated 0-line

    float 
        Current/recalculated 1-line
    
    float
        lower bound precalled purity value
    
    float
        upper bound precalled purity value

    pd.DataFrame
        all samples within the lower and upper bound of the precalled purity value
    """
    # defaults to selecting all precalled purity values
    purity_range_lower = 0
    purity_range_upper = 100

    if step_size is None:
        step_size = 0.05

    data_df = data.df
    data_sample = data_df.loc[data_id] # gets a specific row of sample data
    current_purity_value = data_sample[purity_col] # gets the current purity value
    
    
    if manual_input_source == "Select All":
        
        # displays all precalled purity values!!
        line_0 = 0
        line_1 = 50 

        slider_value = [line_0, line_1]
    else:
        
        # if manual_input_source == "Use slider":
            # slide_value 
            # line_0 = slider_value[0]
            # line_1 = slider_value[1]

        # gets the lower and upper range of purity values
        # accounts for the floating point precision issues
        purity_range_lower = current_purity_value*100 - slider_value  #(line_1 - line_0)
        purity_range_upper = current_purity_value*100 + slider_value #(line_1 - line_0)

        # checks to make sure you dont get negative purity values
        if purity_range_lower < 0:
            purity_range_lower = 0

        # checks to make sure purity value isn't greater than 100%
        if purity_range_upper > 100:
            purity_range_upper = 100
        
    # purity_values = data_df[purity_col] # gets the purity values
    # purity_values_in_range = purity_values.where((purity_values*100 >= purity_range_lower) & (purity_values*100 <= purity_range_upper))
    # purity_value_idx = purity_range_lower
    # precalled_sample_values = data_df.loc[purity_values_in_range]
    
    precalled_sample_values_df = data_df.loc[(data_df[purity_col]*100 >= purity_range_lower) & (data_df[purity_col]*100 <= purity_range_upper)]
    print("precalled sample values dataframe")
    print(precalled_sample_values_df)
    # your parameters have to match the inputs in the same order
    # your returns have to match the outputs in the same order

    return [
        slider_value,
        current_purity_value,
        # purity_range_lower,
        # purity_range_upper,
        # precalled_sample_values
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
        step_size = 5
    return [
            html.Div(
                [
                    # gives user option of filtering out specific purity values based on slider input
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
                                    id="precalled-selection-type-radioitems",
                                ),
                            ])
                        ]),
                    # displays the current purity value
                    dbc.Row([
                        html.Div(
                            [
                                # might need to make a column or do inline component
                                dbc.Label("Current Purity Value: "),
                                html.Label(children="", id="current-purity-value"), # initialize label to empty string
                            ])
                        ]),

                    dbc.Row(
                        [
                            # creating a horizontal slider for selecting a range of purity values
                            # try to find parameter that snaps to whole integer, multiples of 5
                            dcc.Slider(
                                id='custom-precalled-slider', 
                                min=0.0, 
                                max=50.0,
                                step=step_size,
                                value= 5, # initial value on slider 
                                marks={i: f'{i}' for i in range(0, 50, 5)}, 
                                tooltip={"placement": "right", "always_visible": True}
                            )
                        ], 
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

def gen_absolute_precalled_custom_solution_component(step_size=None):
    """
    Generates an AppComponent defining the interactive elements for setting a manual purity/ploidy solution using a copy number profile

    Returns
    =======
    AnnoMate.AppComponent
        AppComponent defining the interactive elements for setting a manual purity/ploidy solution using a copy number profile
    """
    
    return AppComponent(
        'Precalled Purity', # table title
        layout= gen_precalled_absolute_custom_solution_layout(step_size=step_size), # html layout of app component
        new_data_callback=gen_custom_precalled_absolute_component,
        internal_callback=gen_custom_precalled_absolute_component,
        callback_input=[
            Input('custom-precalled-slider', 'value'),
            Input('precalled-selection-type-radioitems', 'value')
            # make sure the order and value for callback_input Input matches the inputs
            # for the gen_custom_precalled_absolute_component
            # i.e. slider input, radio item input
        ],
        callback_output=[
            Output('custom-precalled-slider', 'value'),
            Output('current-purity-value', 'children')


        # slider_value,
        # purity_range_lower,
        # purity_range_upper,
        # precalled_sample_values

            # make sure the order and value for callback_output Output matches the outputs
            # for the gen_custom_precalled_absolute_component
            # i.e. slider input, dataframe of the precalled purity values
        ],
    )
