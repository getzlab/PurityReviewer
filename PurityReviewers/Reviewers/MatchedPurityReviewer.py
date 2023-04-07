from JupyterReviewer.Data import Data, DataAnnotation
from JupyterReviewer.ReviewDataApp import ReviewDataApp, AppComponent
from JupyterReviewer.DataTypes.GenericData import GenericData

import pandas as pd
import numpy as np
import time
import warnings

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
from dash import Dash, dash_table
import dash_bootstrap_components as dbc

from JupyterReviewer.ReviewerTemplate import ReviewerTemplate

from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
import os
import pickle
from typing import Union, List, Dict
import sys

from PurityReviewers.AppComponents.AbsoluteSolutionsReportComponent import gen_absolute_solutions_report_component
from PurityReviewers.AppComponents.AbsoluteCustomSolutionComponent import gen_absolute_custom_solution_component
from PurityReviewers.AppComponents.utils import gen_cnp_figure, gen_mut_figure, parse_absolute_soln, validate_purity, validate_ploidy, csize

from JupyterReviewer.AppComponents.DataTableComponents import gen_annotated_data_info_table_component
    

class MatchedPurityReviewer(ReviewerTemplate):

    def gen_data(self,
                 description: str,
                 df: pd.DataFrame,
                 index: List,
                 # preprocess_data_dir: str,
                 # acs_col,
                 # maf_col,
                 # rdata_fn_col,
                 # reload_cnp_figs=False,
                 # reload_mut_figs=False,
                 # mut_fig_hover_data=[],
                 annot_df: pd.DataFrame = None,
                 annot_col_config_dict: Dict = None,
                 history_df: pd.DataFrame = None) -> GenericData:
        """
        Parameters
        ==========
        preprocess_data_dir: str or path
            Path to directory to store premade plots
        acs_col: str
            Column name in df with path to seg file produced by alleliccapseg
        maf_col: str
            Column name in df with path to mutation validator validated maf
        rdata_fn_col: str
            Column name in df with path to rdata produced by Absolute
        reload_mut_figs: bool
            Whether to regenerate the mutation figures again
        reload_cnp_figs: bool
            Whether to regenerate the copy number plot
        mut_fig_hover_data: List
            List of column names in the maf file (from maf_col in df) to display on hover in
            mutation figures

        """
        pandas2ri.activate()
        # if not os.path.exists(preprocess_data_dir):
        #     os.mkdir(preprocess_data_dir)
        
        # preprocessing
        # 1. download rdata 
#         rdata_dir = f'{preprocess_data_dir}/rdata_to_tsv'
#         if not os.path.exists(rdata_dir):
#             print(f'Converting ABSOLUTE rdata into tsv files in {rdata_dir}')
#             os.mkdir(rdata_dir)
#             df[f'{rdata_fn_col}_as_tsv'] = ''
#             for i, r in df.iterrows():
#                 output_fn = f'{rdata_dir}/{i}.rdata.tsv'
#                 df.loc[i, f'local_absolute_rdata_as_tsv'] = output_fn
#                 try:
#                     parse_absolute_soln(df.loc[i, rdata_fn_col]).to_csv(output_fn, sep='\t')
#                     df.loc[i, f'local_absolute_rdata_as_tsv_downloaded'] = True
#                 except Exception as e:
#                     print(sys.exc_info()[2])
#                     warnings.warn(f'Failed to parse {df.loc[i, rdata_fn_col]}. No tsv available. Original error: {e}')
#                     df.loc[i, f'local_absolute_rdata_as_tsv_downloaded'] = False
#         else:
#             print(f'rdata tsv directory already exists: {rdata_dir}')
#             df[f'{rdata_fn_col}_as_tsv'] = ''
            
#             # existing_files = os.listdir(rdata_dir)
#             # print(existing_files)
#             for i, r in df.iterrows():
#                 output_fn = f'{rdata_dir}/{i}.rdata.tsv'
#                 # if output_fn not in existing_files:
#                 #     raise ValueError(f'Missing file {output_fn}')
                
#                 df.loc[i, f'{rdata_fn_col}_as_tsv'] = output_fn

        # 2. Process cnp figures
        # cnp_figs_dir = f'{preprocess_data_dir}/cnp_figs'
        # if not os.path.exists(cnp_figs_dir):
        #     os.mkdir(cnp_figs_dir)
        #     reload_cnp_figs = True
        # else:
        #     print(f'cnp figs directory already exists: {cnp_figs_dir}')
            
#         if reload_cnp_figs:
#             print('Reloading cnp figs')
#             for i, r in df.iterrows():
#                 output_fn = f'{cnp_figs_dir}/{i}.cnp_fig.pkl'
#                 fig = gen_cnp_figure(df.loc[i, acs_col])
#                 pickle.dump(fig, open(output_fn, "wb"))
#                 df.loc[i, f'cnp_figs_pkl'] = output_fn
                
#         mut_figs_dir = f'{preprocess_data_dir}/mut_figs'
#         if not os.path.exists(mut_figs_dir):
#             os.mkdir(mut_figs_dir)
#             reload_mut_figs = True
#         else:
#             print(f'mut figs directory already exists: {mut_figs_dir}')
      
#         if reload_mut_figs:
#             print('Reloading mut figs')
#             for i, r in df.iterrows():
#                 output_fn = f'{mut_figs_dir}/{i}.cnp_fig.pkl'
#                 fig = gen_mut_figure(df.loc[i, maf_col], hover_data=mut_fig_hover_data)
#                 pickle.dump(fig, open(output_fn, "wb"))
#                 df.loc[i, f'mut_figs_pkl'] = output_fn

        return GenericData(index=index,
                           description=description,
                           df=df,
                           annot_df=annot_df,
                           annot_col_config_dict=annot_col_config_dict,
                           history_df=history_df)

    def gen_review_app(self,
                       sample_info_cols: List[str],
                       acs_col,
                       maf_col,
                       rdata_fn_col,
                       mut_fig_hover_data=[],
                       csize=csize,
                       custom_parse_absolute_soln=None,
                       step_size=None
                       ) -> ReviewDataApp:
        """
        Parameters
        ==========
        sample_info_cols: list[str]
            List of columns in data
        acs_col: str
            Column name in data with path to seg file from alleliccapseg
        maf_col: str
            Column name in data with path to maf file (mutation validator validated maf)
        rdata_fn_col: str
            Column name in data with LOCAL path to maf file. Should be predownloaded at set_review_data()
        cnp_fig_pkl_fn_col: str
            Column nae in data with LOCAL path to premade pickled copy number figures
        mut_fig_pkl_fn_col: str
            Column nae in data with LOCAL path to premade pickled mutation figures
        """

        app = ReviewDataApp()
        
        app.add_component(
            gen_absolute_solutions_report_component(),
            rdata_fn_col=rdata_fn_col,
            acs_col=acs_col, 
            maf_col=maf_col,
            mut_fig_hover_data=mut_fig_hover_data,
            csize=csize,
            custom_parse_absolute_soln=custom_parse_absolute_soln,
        )

        app.add_component(
            gen_annotated_data_info_table_component(), 
            cols=sample_info_cols, 
            data_attribute='df',
            link_display_name=None
        )

        app.add_component(
            gen_absolute_custom_solution_component(step_size=step_size),
            # cnp_fig_pkl_fn_col='cnp_figs_pkl'
            acs_col=acs_col,
            step_size=step_size,
            csize=csize,
        )

        return app
    
    def set_default_autofill(self):
        self.add_autofill('Pick current ABSOLUTE solution', State('absolute-purity', 'children'), 'Purity')
        self.add_autofill('Pick current ABSOLUTE solution', State('absolute-ploidy', 'children'), 'Ploidy')
        self.add_autofill('Pick current ABSOLUTE solution', fill_value='Absolute', annot_name='Method')
        self.add_autofill(
            'Pick current ABSOLUTE solution',
            fill_value=State('absolute-solution-idx', 'children'),
            annot_name='Absolute_solution_idx'
        )
        self.add_autofill('Use current manual solution', State('custom-cnp-graph-purity', 'value'), 'Purity')
        self.add_autofill('Use current manual solution', State('custom-cnp-graph-ploidy', 'value'), 'Ploidy')
        self.add_autofill('Use current manual solution', fill_value='Manual', annot_name='Method')
        

    def set_default_review_data_annotations(self):
        self.add_review_data_annotation(
            annot_name='Purity',
            review_data_annotation=DataAnnotation('float', validate_input=validate_purity))
        self.add_review_data_annotation(
            annot_name='Ploidy',
            review_data_annotation=DataAnnotation('float', validate_input=validate_ploidy))
        
        self.add_review_data_annotation(
            annot_name='Method', 
            review_data_annotation=DataAnnotation(
                annot_value_type='string', 
                options=['Absolute', 'Manual']))

        self.add_review_data_annotation(
            annot_name='Absolute_solution_idx',
            review_data_annotation=DataAnnotation(annot_value_type='int'))
        
        self.add_review_data_annotation(
            annot_name='Notes', 
            review_data_annotation=DataAnnotation(
                annot_value_type='string')
        )

    def set_default_review_data_annotations_app_display(self):
        self.add_review_data_annotations_app_display(annot_name='Purity', app_display_type='number')
        self.add_review_data_annotations_app_display(annot_name='Ploidy', app_display_type='number')
        self.add_review_data_annotations_app_display(annot_name='Method', app_display_type='select')
        self.add_review_data_annotations_app_display(annot_name='Absolute_solution_idx', app_display_type='number')
        self.add_review_data_annotations_app_display(annot_name='Notes', app_display_type='textarea')

        