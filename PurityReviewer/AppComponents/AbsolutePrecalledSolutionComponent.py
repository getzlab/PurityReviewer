"""
Displays a allelic copy ratio profile with the option to set the 0 and 1 line (via slider or input value) corresponding to an integer assignment to copy number peaks. Automatically calculates the purity given the corresponding solution.
"""
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from cnv_suite.visualize import plot_acr_interactive
from PurityReviewer.AppComponents.utils import gen_cnp_figure, gen_mut_figure, CSIZE_DEFAULT, parse_absolute_soln, calculate_multiplicity

from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
import os
import pickle
from typing import Union, List, Dict
import sys
from cnv_suite import calc_cn_levels
import pandas as pd
import numpy as np

MANUAL_INPUT_SOURCE = ["Use slider", "Select All"]
absolute_rdata_cols = ['alpha', 'tau', 'tau_hat', '0_line', '1_line',
                       'sigma_H', 
                       'theta_Q', 
                       'lambda',  
                       'SCNA_likelihood', 
                       'Kar_likelihood', 
                       'SSNVs_likelihood']
SELECTED_ROW_INDEX = 2
CURRENT_ABSOLUTE_SOLUTION_IDX = -1
ABSOLUTE_PURITY_VALUE_COLUMN_NAME = 'alpha'

# Add a custom component: below, I add a component that allows you to manually set the range of purity values
# NOTE: the purity calculated is tau, NOT tau_g. 
# Use the called purity and tau as inputs to absolute_segforcecall to get tau_g (tau_hat)
def gen_custom_precalled_absolute_component(
    data: GenericData,
    data_id,
    slider_value,  # dash app parameters come first
    selected_row_array, # dash app parameters come first
    manual_input_source,
    rdata_fn_col,
    acs_col, 
    mut_fig_hover_data,
    purity_col='PCA__ABSOLUTE__Cancer_DNA_fraction',
    step_size=None,
    custom_parse_absolute_soln=None
):
    """
    Callback function that updates copy number plots when a new sample/pair is selected 
    or a new purity value range is indicated via slider or select all radio item
    
    Parameters
    ==========
    data: GenericData
        Data object storing the relevant data to review
        
    data_id: str
        Name of the object being reviewed
        
    slider_value: List[float]
        List of length 2 where the 0-index entry is the current value of the 0-line, and the 1-index entry is the current value of the 1-line

    selected_row_array: List
        List of length 1 containing the currently selected ABSOLUTE solution from the ABSOLUTE solution table
        
    manual_input_source: str ["Use slider", "Select All"]
        Which mode to set get the purity solution and replot the copy number profile with the corresponding comb

    rdata_fn_col: str
        Column in data.df corresponding to the LOCAL file path of the ABSOLUTE rdata

    acs_col: str
        Column name in data with path to seg file from alleliccapseg or other tsv with allelic copy ratio measurements.
        
    mut_fig_hover_data: List[str]
        List of column names to add to plotly hover_data in mutation figure
        
    purity_col: str
        Column name in data with path to precalled purity values

    step_size: int, default=5
        Degree of precision for the slider bar
    
    custom_parse_absolute_soln: function
        Custom absolute parser function (rdata_path -> data_df)

    Returns
    =======
    List[float]
        List of length 2 where the 0-index entry is the current/recalculated value of the 0-line, and the 1-index entry is the current value of the 1-line. Used to update the slider bar

    float
        Current/recalculated purity value
    
    List[int]
        An array of length 1 indicating the selected row in the ABSOLUTE solution table

    List[Dict]
        ABSOLUTE solutions table converted into a list of records to be displayed in the dashboard 
        with only purity values within a predefined range of the current purity value
    
    plotly.Figure
        Copy number profile plot with the "comb" plotted corresponding to the currently selected solution

    plotly.Figure
        Mutation profile plot with the "comb" plotted corresponding to the currently selected solution

    float
        Purity corresponding to the currently selected solution. For new data, it is set to the first solution

    float
        Ploidy corresponding to the currently selected solution. For new data, it is set to the first solution

    int
        Index of the current ABSOLUTE solution. For new data, it is set to the first solution 0
    """    
    # defaults to selecting all precalled purity values
    purity_range_lower = 0
    purity_range_upper = 100

    if step_size is None:
        step_size = 5

    data_df = data.df.copy() # do not want to modify the original pairs dataframe 
    data_sample = data_df.loc[data_id] # gets a specific row of paired sample data
    current_purity_value = data_sample[purity_col] # gets the current tumor sample purity value

    # checks if you want to limit the number of samples displayed of purity values
    if manual_input_source == "Use slider":
    
        # gets the lower and upper range of purity values
        # accounts for the floating point precision issues
        purity_range_lower = current_purity_value*100 - slider_value  
        purity_range_upper = current_purity_value*100 + slider_value

        # checks to make sure you dont get negative purity values
        if purity_range_lower < 0:
            purity_range_lower = 0

        # checks to make sure purity value isn't greater than 100%
        if purity_range_upper > 100:
            purity_range_upper = 100
        
    # from absolute solutions report component gen_absolute_solutions_report_new_data
    parse_absolute_soln_func = custom_parse_absolute_soln if custom_parse_absolute_soln is not None else parse_absolute_soln
    try:
        absolute_rdata_df,maf,maf_annot_list = parse_absolute_soln_func(data_sample[rdata_fn_col])
    except Exception as e:
        print(e)
        print("excepted error and initialized empty dataframes")
        absolute_rdata_df,maf,maf_annot_list = pd.DataFrame()

    absolute_rdata_df = absolute_rdata_df.round(2)
    
    cnp_fig = gen_cnp_figure(data_sample[acs_col], csize=CSIZE_DEFAULT)

    # add 1 and 0 lines
    cnp_fig_with_lines = go.Figure(cnp_fig)
    
    purity = 0
    ploidy = 0

    absolute_rdata_within_range_df = absolute_rdata_df.loc[(absolute_rdata_df[ABSOLUTE_PURITY_VALUE_COLUMN_NAME]*100 <= purity_range_upper) 
                                                         & (absolute_rdata_df[ABSOLUTE_PURITY_VALUE_COLUMN_NAME]*100 >= purity_range_lower)]
    
    if absolute_rdata_within_range_df.shape[0] > 0:
        solution_data = absolute_rdata_within_range_df.iloc[selected_row_array[0]]
        maf_soln = pd.concat([maf, maf_annot_list[selected_row_array[0]]],axis=1)
        maf_soln = maf_soln[maf_soln['Variant_Type']=='SNP']

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
        
        purity = solution_data['alpha']
        ploidy = solution_data['tau_hat']

        maf_soln['multiplicity'] = calculate_multiplicity(maf_soln,purity)
        mut_fig = gen_mut_figure(maf_soln, hover_data=mut_fig_hover_data, csize=CSIZE_DEFAULT)
        mut_fig_with_lines = go.Figure(mut_fig)
        for yval in [1,2]:
            mut_fig_with_lines.add_hline(y=yval,
                                    line_dash="dash",
                                    line_color='black',
                                    line_width=1)
     
    return [
        slider_value,
        current_purity_value,
        selected_row_array,
        absolute_rdata_within_range_df.to_dict('records'),
        cnp_fig_with_lines, 
        mut_fig_with_lines,
        purity,
        ploidy, 
        0
    ]

def gen_precalled_solutions_report_internal(
    data: GenericData,
    data_id,
    slider_value,  # dash app parameters come first
    selected_row_array, # dash app parameters come first
    manual_input_source,
    rdata_fn_col,
    acs_col, 
    mut_fig_hover_data,
    purity_col='PCA__ABSOLUTE__Cancer_DNA_fraction',
    step_size=None,
    custom_parse_absolute_soln=None
):
    """
    Callback function that updates copy number plots when the selected ABSOLUTE solution changes
    
    Parameters
    ==========
    data: GenericData
        Data object storing the relevant data to review
        
    data_id: str
        Name of the object being reviewed
        
    selected_row_array: List
        List of length 1 containing the currently selected ABSOLUTE solution from the ABSOLUTE solution table
        
    rdata_fn_col: str
        Column in data.df corresponding to the LOCAL file path of the ABSOLUTE rdata
        
    acs_col: str
        Column name in data with path to seg file from alleliccapseg or other tsv with allelic copy ratio measurements.
        
    mut_fig_hover_data: List[str]
        List of column names to add to plotly hover_data in mutation figure

    purity_col: str
        Column name in data with path to precalled purity values

    step_size: int, default=5
        Degree of precision for the slider bar
            
    custom_parse_absolute_soln: function
        Custom absolute parser function (rdata_path -> data_df)

    Returns
    =======
    List[float]
        List of length 2 where the 0-index entry is the current/recalculated value of the 0-line, and the 1-index entry is the current value of the 1-line. Used to update the slider bar

    float
        Current/recalculated purity value
    
    List[int]
        An array of length 1 indicating the selected row in the ABSOLUTE solution table

    List[Dict]
        ABSOLUTE solutions table converted into a list of records to be displayed in the dashboard 
        with only purity values within a predefined range of the current purity value
    
    plotly.Figure
        Copy number profile plot with the "comb" plotted corresponding to the currently selected solution

    plotly.Figure
        Mutation profile plot with the "comb" plotted corresponding to the currently selected solution

    float
        Purity corresponding to the currently selected solution. For new data, it is set to the first solution

    float
        Ploidy corresponding to the currently selected solution. For new data, it is set to the first solution

    int
        Index of the current ABSOLUTE solution. For new data, it is set to the first solution 0
    """
    
    output_data = gen_custom_precalled_absolute_component(
        data,
        data_id,
        slider_value,  # dash app parameters come first
        selected_row_array, # dash app parameters come first
        manual_input_source,
        rdata_fn_col,
        acs_col, 
        mut_fig_hover_data=mut_fig_hover_data,
        step_size=step_size,
        purity_col=purity_col,
        custom_parse_absolute_soln=custom_parse_absolute_soln
    )

    output_data[SELECTED_ROW_INDEX] = selected_row_array
    output_data[CURRENT_ABSOLUTE_SOLUTION_IDX] = selected_row_array[0] + 1 # 1 based indexing for copy number profile

    return output_data
    
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
            # displays the interactive component to filter the samples displays based on their purity values
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
                    ),
                    # displays the absolute solutions report
                html.Div(
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
            ],
                  
        )
    ]

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
        internal_callback=gen_precalled_solutions_report_internal,
        callback_input=[
            Input('custom-precalled-slider', 'value'),
            Input('absolute-rdata-select-table', 'selected_rows'),
            Input('precalled-selection-type-radioitems', 'value'),
        ],
        callback_output=[
            Output('custom-precalled-slider', 'value'),
            Output('current-purity-value', 'children'),
            Output('absolute-rdata-select-table', 'selected_rows'),
            Output('absolute-rdata-select-table', 'data'),
            Output('cnp-graph', 'figure'),
            Output('mut-graph', 'figure'),
            Output('absolute-purity', 'children'),
            Output('absolute-ploidy', 'children'),
            Output('absolute-solution-idx', 'children'),
        ],
    )