"""
Given a precalled purity sample value displays the absolute solutions within a specific 
range of the current precalled purity value. The default range is 5% above and below 
the current precalled purity value. Can also choose to see all absolute solutions.
"""
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from cnv_suite.visualize import plot_acr_interactive
from PurityReviewer.AppComponents.utils import gen_cnp_figure, gen_mut_figure, CSIZE_DEFAULT, parse_absolute_soln, calculate_multiplicity, mut_af_plot, multiplicity_plot
#gen_allele_fraction_figure

from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
import os
import pickle
from typing import Union, List, Dict
import sys
from cnv_suite import calc_cn_levels
import pandas as pd
import numpy as np

PRECALLED_SLIDER_VALUES = ["Use slider", "Select All"]
absolute_rdata_cols = ['alpha', 'tau', 'tau_hat', '0_line', '1_line',
                       'sigma_H', 
                       'theta_Q', 
                       'lambda',  
                       'SCNA_likelihood', 
                       'Kar_likelihood', 
                       'SSNVs_likelihood']
PRECALLED_SLIDER_VALUE_INDEX = 0
SELECTED_ROW_INDEX = 2
CURRENT_ABSOLUTE_SOLUTION_IDX = -1
ABSOLUTE_PURITY_VALUE_COLUMN_NAME = 'alpha'
PRECALLED_PURITY_COLUMN_NAME = 'precalled_purity_values'

# NOTE: the purity calculated is tau, NOT tau_g. 
# Use the called purity and tau as inputs to absolute_segforcecall to get tau_g (tau_hat)
def gen_absolute_solutions_report_range_of_precalled_component( 
    data: GenericData,
    data_id,
    slider_value,  # dash app parameters come first
    selected_row_array, # dash app parameters come first
    precalled_radio_button_value,
    rdata_fn_col,
    acs_col, 
    mut_fig_hover_data,
    custom_parse_absolute_soln=None,
    internal_callback=False
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
        
    slider_value: int
        the value for the range of purity values to view centered around the precalled purity value

    selected_row_array: List[int]
        List of length 1 containing the currently selected ABSOLUTE solution from the ABSOLUTE solution table
        
    precalled_radio_button_value: str ["Use slider", "Select All"]
        Which mode the user wants to use for filtering the absolute solutions that are viewed, based 
        on the precalled purity value; either looking at all the absolute solutions or looking 
        absolute solutions within specific range of precalled purity value  

    rdata_fn_col: str
        Column in data.df corresponding to the LOCAL file path of the ABSOLUTE rdata

    acs_col: str
        Column name in data with path to seg file from alleliccapseg or other tsv with allelic copy ratio measurements.
        
    mut_fig_hover_data: List[str]
        List of column names to add to plotly hover_data in mutation figure
    
    custom_parse_absolute_soln: function
        Custom absolute parser function (rdata_path -> data_df)
    
    internal_callback: boolean, default is False
        whether the function is being called from the new_callback or internal_callback

    Returns
    =======
    int
        the value for the range of purity values to view centered around the precalled purity value

    float
        precalled purity value
    
    List[int]
        An array of length 1 indicating the selected row in the ABSOLUTE solution table

    List[Dict]
        ABSOLUTE solutions table converted into a list of records to be displayed in the dashboard 
        with only purity values within a predefined range of the precalled purity value
    
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
    print("inside the newdata callback function")


    # checks if you are loading in the data for the first time (new data callback)
    if not internal_callback:
        slider_value = 5

    # defaults to selecting all precalled purity values
    purity_range_lower = 0
    purity_range_upper = 100

    data_df = data.df.copy() # do not want to modify the original pairs dataframe 
    pairs_data_row = data_df.loc[data_id] # gets a specific row of paired sample data
    precalled_purity_value = pairs_data_row[PRECALLED_PURITY_COLUMN_NAME] # gets the precalled purity value from tumor sample

    # checks if you want to display less absolute solutions based on its purity value 
    if precalled_radio_button_value == "Use slider":
    
        # gets the lower and upper range of purity values
        purity_range_lower = precalled_purity_value*100 - slider_value  
        purity_range_upper = precalled_purity_value*100 + slider_value

        # checks to make sure you dont get negative purity values
        if purity_range_lower < 0:
            purity_range_lower = 0

        # checks to make sure purity value isn't greater than 100%
        if purity_range_upper > 100:
            purity_range_upper = 100
        
    # from absolute solutions report component gen_absolute_solutions_report_new_data
    parse_absolute_soln_func = custom_parse_absolute_soln if custom_parse_absolute_soln is not None else parse_absolute_soln
    try:
        absolute_rdata_df,maf,maf_annot_list = parse_absolute_soln_func(pairs_data_row[rdata_fn_col])
    except Exception as e:
        print(e)
        print("excepted error and initialized empty dataframes")
        absolute_rdata_df,maf,maf_annot_list = pd.DataFrame()

    absolute_rdata_df = absolute_rdata_df.round(2)
    
    cnp_fig = gen_cnp_figure(pairs_data_row[acs_col], csize=CSIZE_DEFAULT)

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

        # maf variable contains all the mutation data, the maf_df

        # if you want to find the mutation maf data use the maf_soln dataframe !!! and the alt_count and ref_count columns naming convention presumably will remain the same
        # check column names of maf_soln because in the calculate_multiplicity function the alt_count column is accessed with the key 'alt'
        maf_soln['multiplicity'] = calculate_multiplicity(maf_soln, purity)

        mut_fig = gen_mut_figure(maf_soln, hover_data=mut_fig_hover_data, csize=CSIZE_DEFAULT)
        mut_fig_with_lines = go.Figure(mut_fig)
        # ssnv_colors = ["dodgerblue", "darkgrey", "seagreen3"]
        # draw_indvividual_lines = True
        
        # allele_fraction_fig = mut_af_plot(maf_soln, ssnv_colors, draw_indvividual_lines) 
        # {"af_post_pr": allele_frac_post_probability, "grid_mat": grid_mat}
        # allele_fraction_fig, af_probabilities_grid_dict = mut_af_plot(maf_soln, ssnv_colors, draw_indvividual_lines) 
        allele_fraction_fig, _, debugging_value = mut_af_plot(maf_soln, 
                                            #  ssnv_colors, draw_indvividual_lines
                                             ) 

        # # DEBUGGING !!!
        # debugging_value = f"shape: {debugging.shape} values: "
        # for val in debugging[:5]:
        #     debugging_value = debugging_value + "" + str(val) + ", " 
        # # END OF DEBUGGING!!


        # allele_frac_posterior_probability = af_probabilities_grid_dict['af_post_pr']
        # grid_mat = af_probabilities_grid_dict['grid_mat'] # NEED TO FIGURE OUT WHAT THIS DOES
        
        mode_color = "grey"
        draw_indv = True
        # ssnv_multiplicity_fig = multiplicity_plot(maf_soln, allele_frac_posterior_probability, grid_mat, mode_color, draw_indv)
        # allele_fraction_fig = gen_allele_fraction_figure(maf_soln) 
        allele_fraction_lines = go.Figure(allele_fraction_fig)
        # ssnv_multiplicity_lines = go.Figure()
        # ssnv_multiplicity_lines = go.Figure(ssnv_multiplicity_fig)

        for yval in [1,2]:
            mut_fig_with_lines.add_hline(y=yval,
                                    line_dash="dash",
                                    line_color='black',
                                    line_width=1)
    
    return [
        5, # default range for the slider
        precalled_purity_value,
        [0], # default selected row is first row in table
        absolute_rdata_within_range_df.to_dict('records'),
        cnp_fig_with_lines, 
        mut_fig_with_lines,
        allele_fraction_fig,
        debugging_value,
        # allele_fraction_lines,
        # ssnv_multiplicity_lines,
        purity,
        ploidy, 
        1 # defaults to having the 1st copy number profile 
    ]

def gen_absolute_precalled_solutions_report_internal(
    data: GenericData,
    data_id,
    slider_value,  # dash app parameters come first
    selected_row_array, # dash app parameters come first
    precalled_radio_button_value,
    rdata_fn_col,
    acs_col, 
    mut_fig_hover_data,
    custom_parse_absolute_soln=None,
    internal_callback=False
):
    """
    Callback function that updates copy number plots when the selected ABSOLUTE solution changes
    
    Parameters
    ==========
    data: GenericData
        Data object storing the relevant data to review
        
    data_id: str
        Name of the object being reviewed

    slider_value: int
        the value for the range of purity values to view centered around the precalled purity value
        
    selected_row_array: List
        List of length 1 containing the currently selected ABSOLUTE solution from the ABSOLUTE solution table
        
    precalled_radio_button_value: str ["Use slider", "Select All"]
        Which mode the user wants to use for filtering the absolute solutions that are viewed, based 
        on the precalled purity value; either looking at all the absolute solutions or looking 
        absolute solutions within specific range of precalled purity value 

    rdata_fn_col: str
        Column in data.df corresponding to the LOCAL file path of the ABSOLUTE rdata
        
    acs_col: str
        Column name in data with path to seg file from alleliccapseg or other tsv with allelic copy ratio measurements.
        
    mut_fig_hover_data: List[str]
        List of column names to add to plotly hover_data in mutation figure
            
    custom_parse_absolute_soln: function
        Custom absolute parser function (rdata_path -> data_df)
    
    internal_callback: boolean, default is False
        whether the function is being called from the new_callback or internal_callback

    Returns
    =======
    int
        the value for the range of purity values to view centered around the precalled purity value

    float
        precalled purity value
    
    List[int]
        An array of length 1 indicating the selected row in the ABSOLUTE solution table

    List[Dict]
        ABSOLUTE solutions table converted into a list of records to be displayed in the dashboard 
        with only purity values within a predefined range of the precalled purity value
    
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
    output_data = gen_absolute_solutions_report_range_of_precalled_component(
        data,
        data_id,
        slider_value,# dash app parameters come first
        selected_row_array, 
        precalled_radio_button_value,
        rdata_fn_col,
        acs_col, 
        mut_fig_hover_data,
        custom_parse_absolute_soln=custom_parse_absolute_soln,
        internal_callback=True
    )

    output_data[PRECALLED_SLIDER_VALUE_INDEX] = slider_value
    output_data[SELECTED_ROW_INDEX] = selected_row_array
    output_data[CURRENT_ABSOLUTE_SOLUTION_IDX] = selected_row_array[0] + 1 # 1 based indexing for copy number profile

    return output_data
    
def gen_absolute_precalled_solution_report_layout():
    """
    Generates the layout of the component in the dashboard corresponding to the report of the 
    ABSOLUTE solutions within a user-defined range of the precalled purity value

    Returns
    =======
    dash.html
        a plotly dash layout with a slider to define the range of purity values, Table with 
        selectable rows for the ABSOLUTE solutions, a copy number profile, and mutation profile
    """
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
                                        } for v in PRECALLED_SLIDER_VALUES
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
                                marks={i: f'{i}' for i in range(0, 51, 5)}, 
                                tooltip={"placement": "right", "always_visible": True}
                            )
                        ], 
                    ),

                # DEBUGGING ELEMENTS
                html.Div(
                    children=[
                        dbc.Label(children='Debugging Some Stuff', id='debugging-value-1')
                    ]
                ),

                # END OF DEBUGGING ELEMENTS!!!

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
                        dcc.Graph(id='mut-graph', figure={}),
                        # allele fraction graph
                        dbc.Label("LOCATION OF THE ALLELE FRACTION GRAPH"),
                        dcc.Graph(id='allele-fraction-graph', figure={}),
                        # multiplicity graph
                        dbc.Label("LOCATION OF THE MULTIPLICITY GRAPH"),
                        dcc.Graph(id='ssnv-multiplicity-graph', figure={}),
                    ]
                )
            ],       
        )
    ]

def gen_absolute_precalled_solutions_report_component():
    """
    Generates an AppComponent defining the interactive elements for viewing ABSOLUTE solutions 
    with purity values within a user-defined range of the precalled purity value
    
    Returns
    =======
    AnnoMate.AppComponent
        AppComponent defining the interactive elements for viewing ABSOLUTE solutions
        with purity values within a user-defined range of the precalled purity value
    """
    
    return AppComponent(
        'Precalled Purity', # table title
        layout= gen_absolute_precalled_solution_report_layout(), # html layout of app component
        new_data_callback=gen_absolute_solutions_report_range_of_precalled_component,
        internal_callback=gen_absolute_precalled_solutions_report_internal,
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
            Output('allele-fraction-graph', 'figure'),
            Output('debugging-value-1', 'children'),
            # Output('ssnv-multiplicity-graph', 'figure'),
            Output('absolute-purity', 'children'),
            Output('absolute-ploidy', 'children'),
            Output('absolute-solution-idx', 'children'),
        ],
    )