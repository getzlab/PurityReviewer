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
SELECTED_ROW_INDEX = 4
CURRENT_ABSOLUTE_SOLUTION_IDX = -1

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
    maf_col,
    mut_fig_hover_data,
    purity_col='PCA__ABSOLUTE__Cancer_DNA_fraction',
    step_size=None,
    csize=None,
    custom_parse_absolute_soln=None
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
        step_size = 5

    data_df = data.df.copy() # do not want to modify the original pairs dataframe 
    data_sample = data_df.loc[data_id] # gets a specific row of sample data
    current_purity_value = data_sample[purity_col] # gets the current purity value

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
        
    # FOUND THE ERROR!!
        # there are more sample values in the data_df dataframe then the absolute_rdata_df!!!!
        # i.e. len(purity_values) > len(absolute_rdata_df) 
    purity_values = data_df[purity_col].copy() # gets the purity values
    purity_values = purity_values.reset_index(drop=True) # creates numerical indices and drops original index
    print("length of purity_values: ", len(purity_values))
    # gets indices of samples with purity values within range
    indices_of_purity_values_in_range = np.where((purity_values*100 >= purity_range_lower) & (purity_values*100 <= purity_range_upper))[0] # makes sure you get the array of indices

    # DEBUGGING!!
    print("these are indices that are within range of the current purity value")
    print(indices_of_purity_values_in_range)
    print() 
    # print(data_df[purity_col].iloc[indices_of_purity_values_in_range])

    # from absolute solutions report component gen_absolute_solutions_report_new_data
    parse_absolute_soln_func = custom_parse_absolute_soln if custom_parse_absolute_soln is not None else parse_absolute_soln
    try:
        absolute_rdata_df,maf,maf_annot_list = parse_absolute_soln_func(data_sample[rdata_fn_col])
    except Exception as e:
        print(e)
        print("except error: and initialize empty dataframes")
        absolute_rdata_df,maf,maf_annot_list = pd.DataFrame()

    absolute_rdata_df = absolute_rdata_df.round(2)
    print("this is the length of parse_absolute_soln: ", len(absolute_rdata_df))
    print("these are the indices for absolute_rdata_df: ", list(absolute_rdata_df.index))
    
    cnp_fig = gen_cnp_figure(data_sample[acs_col], csize=CSIZE_DEFAULT)

    # add 1 and 0 lines
    cnp_fig_with_lines = go.Figure(cnp_fig)
    
    purity = 0
    ploidy = 0
    
    if absolute_rdata_df.shape[0] > 0:
        solution_data = absolute_rdata_df.iloc[selected_row_array[0]]
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

        #mut_fig_with_lines.update_yaxes(range=[0, half_1_line * 2])
        
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
    
    # removes the indices that are out of range error, since len(purity_values) > len(absolute_rdata_df)
    indices_within_bound = np.where(indices_of_purity_values_in_range <= np.max(absolute_rdata_df.index), True, False)
    indices_of_purity_values_in_range = indices_of_purity_values_in_range[indices_within_bound]

    absolute_rdata_within_range_df = absolute_rdata_df.iloc[indices_of_purity_values_in_range]

    for col_nm in absolute_rdata_cols:
    # absolute_rdata_within_range_df = data_df.loc[indices_of_purity_values_in_range, absolute_rdata_cols]
        print(f"checking if column -> {col_nm} is in data_df: {col_nm in data_df.columns}")
    
    # DEBUGGING!!!
    # precalled_sample_values_df = data_df.loc[(data_df[purity_col]*100 >= purity_range_lower) & (data_df[purity_col]*100 <= purity_range_upper)]
    # print("precalled sample values dataframe")
    # print(precalled_sample_values_df)

    return [
        slider_value,
        current_purity_value,
        purity_range_lower,
        purity_range_upper,
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
    maf_col,
    mut_fig_hover_data,
    purity_col='PCA__ABSOLUTE__Cancer_DNA_fraction',
    step_size=None,
    csize=None,
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
        
    maf_col: str
        Column name in data with path to maf file (mutation validator validated maf)
        
    mut_fig_hover_data: List[str]
        List of column names to add to plotly hover_data in mutation figure

    csize: dict
            Dictionary with chromosome sizes (see AppComponents.utils.CSIZE_DEFAULT for hg19)
            
    custom_parse_absolute_soln: function
        Custom absolute parser function (rdata_path -> data_df)

    Returns
    =======
    List[Dict]
        ABSOLUTE solutions table converted into a list of records to be displayed in the dashboard

    plotly.Figure
        Copy number profile plot with the "comb" plotted corresponding to the currently selected solution

    plotly.Figure
        Mutation profile plot with the "comb" plotted corresponding to the currently selected solution

    float
        Purity corresponding to the currently selected solution

    float
        Ploidy corresponding to the currently selected solution

    List[int]
        An array of length 1 indicating the selected row in the ABSOLUTE solution table
    int
        Index of the current ABSOLUTE solution
    """
    
    output_data = gen_custom_precalled_absolute_component(
        data,
        data_id,
        slider_value,  # dash app parameters come first
        selected_row_array, # dash app parameters come first
        manual_input_source,
        rdata_fn_col,
        acs_col, 
        maf_col,
        mut_fig_hover_data,
        csize=CSIZE_DEFAULT,
        custom_parse_absolute_soln=custom_parse_absolute_soln
    )

    output_data[SELECTED_ROW_INDEX] = selected_row_array
    output_data[CURRENT_ABSOLUTE_SOLUTION_IDX] = selected_row_array[0] + 1 # 1 based indexing

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

                    # displays the lower bound for the purity values that will be displayed
                    dbc.Row([
                        html.Div(
                            [
                                dbc.Label("Lower bound range of Purity Values to display: "),
                                html.Label(children="", id="lower-bound-of-purity-values"), # initialize label to empty string
                            ])
                        ]),

                    # displays the upper bound for the purity values that will be displayed
                    dbc.Row([
                        html.Div(
                            [
                                dbc.Label("Upper bound range of Purity Values to display: "),
                                html.Label(children="", id="upper-bound-of-purity-values"), # initialize label to empty string
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

# def filter_purity_values(data, purity_val_lower_range, purity_val_upper_range, 
#                          current_purity_val, purity_val_col) -> GenericData:
#     """
#     Only returns the purity values within a specific range

#     Parameters:
#         data: GenericData, 
#         purity_value_lower_range:
#         purity_value_upper_range:
#         purity_val_col: str, name of the column that has the purity value

#     Return:
#         returns a generic data object with only rows that contained purity values 
#     in a specific range
#     """
#     data = data.copy()

#     data_df = data.df
#     annotation_df = data.annot_df
#     purity_val_range = purity_val_upper_range - purity_val_lower_range

#     filter_row_idx = pd.where((data_df[purity_val_col] >= current_purity_val - purity_val_range) and (data_df[purity_val_col] <= current_purity_val + purity_val_range))

#     # only gets the purity values within a specific range of the current purity value
#     data.df = data_df.iloc[filter_row_idx]
#     data.annot_df = annotation_df.iloc[filter_row_idx]

#     return data

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
            Input('absolute-rdata-select-table', 'selected_rows'),
            Input('precalled-selection-type-radioitems', 'value'),
        ],
        callback_output=[
            Output('custom-precalled-slider', 'value'),
            Output('current-purity-value', 'children'),
            Output('lower-bound-of-purity-values', 'children'),
            Output('upper-bound-of-purity-values', 'children'),
            Output('absolute-rdata-select-table', 'selected_rows'),
            Output('absolute-rdata-select-table', 'data'),
            Output('cnp-graph', 'figure'),
            Output('mut-graph', 'figure'),
            Output('absolute-purity', 'children'),
            Output('absolute-ploidy', 'children'),
            Output('absolute-solution-idx', 'children'),
        ],
    )