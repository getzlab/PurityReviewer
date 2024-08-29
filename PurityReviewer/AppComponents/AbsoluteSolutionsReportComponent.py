"""
Displays a allelic copy ratio profile and a table of solutions from ABSOLUTE (Carter, 2014). The rows of the table can be selected, and the copy number profile plot will be updated with the corresponding ABSOLUTE "comb" solution. 
"""
import pandas as pd
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import plotly.graph_objects as go

from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from cnv_suite.visualize import plot_acr_interactive
from PurityReviewer.AppComponents.utils import gen_cnp_figure, gen_mut_figure, CSIZE_DEFAULT, parse_absolute_soln, calculate_multiplicity


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
    """
    Callback function that updates copy number plots when a new sample/pair is selected 
    
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
        Purity corresponding to the currently selected solution. For new data, it is set to the first solution

    float
        Ploidy corresponding to the currently selected solution. For new data, it is set to the first solution

    List[int]
        An array of length 1 indicating the selected row in the ABSOLUTE solution table. For new data, it is set to the first solution [0]
    int
        Index of the current ABSOLUTE solution. For new data, it is set to the first solution 0
    """
    
    data_df = data.df
    r = data_df.loc[data_id]
    
    parse_absolute_soln_func = custom_parse_absolute_soln if custom_parse_absolute_soln is not None else parse_absolute_soln
    try:
        absolute_rdata_df,maf,maf_annot_list = parse_absolute_soln_func(r[rdata_fn_col])
    except Exception as e:
        print(e)
        absolute_rdata_df,maf,maf_annot_list = pd.DataFrame()

    absolute_rdata_df = absolute_rdata_df.round(2)
    
    cnp_fig = gen_cnp_figure(r[acs_col], csize=CSIZE_DEFAULT)



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



#
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
    """
    Generates the layout of the ABSOLUTE solutions report component in the dashboard

    Returns
    =======
    dash.html
        a plotly dash layout with a Table with selectable rows for the ABSOLUTE solutions, a copy number profile, and mutation profile
    """
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
    """
    Generates an AppComponent defining the interactive elements for viewing ABSOLUTE solutions

    Returns
    =======
    AnnoMate.AppComponent
        AppComponent defining the interactive elements for viewing ABSOLUTE solutions
    """
    
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
    
    
    