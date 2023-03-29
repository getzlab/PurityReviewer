from JupyterReviewer.Data import Data, DataAnnotation
from JupyterReviewer.ReviewDataApp import ReviewDataApp, AppComponent
from JupyterReviewer.DataTypes.GenericData import GenericData

import pandas as pd
import numpy as np

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

from JupyterReviewer.ReviewerTemplate import ReviewerTemplate
from cnv_suite.visualize import plot_acr_interactive

from rpy2.robjects import r, pandas2ri
import os
import pickle
from typing import List, Dict

from PurityReviewers.AppComponents.AbsoluteCustomSolutionComponent import gen_absolute_custom_solution_component
from PurityReviewers.AppComponents.utils import gen_cnp_figure, gen_mut_figure, parse_absolute_soln, validate_purity, validate_ploidy

from JupyterReviewer.AppComponents.DataTableComponents import gen_annotated_data_info_table_component


class ManualPurityReviewer(ReviewerTemplate):

    def gen_data(self,
                 description: str,
                 df: pd.DataFrame,
                 index: List,
                 preprocess_data_dir: str,
                 acs_col,
                 maf_col,
                 rdata_fn_col,
                 reload_cnp_figs=False,
                 reload_mut_figs=False,
                 mut_fig_hover_data=[],
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

        Returns
        =======
        GenericData
            A data object
        """
        pandas2ri.activate()
        if not os.path.exists(preprocess_data_dir):
            os.mkdir(preprocess_data_dir)

        # 2. Process cnp figures
        cnp_figs_dir = f'{preprocess_data_dir}/cnp_figs'
        if not os.path.exists(cnp_figs_dir):
            os.mkdir(cnp_figs_dir)
            reload_cnp_figs = True
        else:
            print(f'cnp figs directory already exists: {cnp_figs_dir}')

        if reload_cnp_figs:
            print('Reloading cnp figs')
            for i, r in df.iterrows():
                output_fn = f'{cnp_figs_dir}/{i}.cnp_fig.pkl'
                fig = gen_cnp_figure(df.loc[i, acs_col])
                pickle.dump(fig, open(output_fn, "wb"))
                df.loc[i, f'cnp_figs_pkl'] = output_fn

        mut_figs_dir = f'{preprocess_data_dir}/mut_figs'
        if not os.path.exists(mut_figs_dir):
            os.mkdir(mut_figs_dir)
            reload_mut_figs = True
        else:
            print(f'mut figs directory already exists: {mut_figs_dir}')


        return GenericData(index=index,
                           description=description,
                           df=df,
                           annot_df=annot_df,
                           annot_col_config_dict=annot_col_config_dict,
                           history_df=history_df)

    def gen_review_app(self,
                       sample_info_cols,
                       acs_col,
                       maf_col,
                       rdata_tsv_fn='local_absolute_rdata_as_tsv',
                       cnp_fig_pkl_fn_col='cnp_figs_pkl',
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
        rdata_tsv_fn: str
            Column name in data with LOCAL path to maf file. Should be predownloaded at set_review_data()
        cnp_fig_pkl_fn_col: str
            Column nae in data with LOCAL path to premade pickled copy number figures
        mut_fig_pkl_fn_col: str
            Column nae in data with LOCAL path to premade pickled mutation figures
        """

        app = ReviewDataApp()
        
        app.add_component(
            gen_annotated_data_info_table_component(), 
            cols=sample_info_cols, 
            data_attribute='df'
        )

        app.add_component(
            gen_absolute_custom_solution_component(step_size=step_size),
            cnp_fig_pkl_fn_col='cnp_figs_pkl',
            step_size=step_size
        )
        
        return app

    def set_default_autofill(self):
        self.add_autofill('cnp-plot-button', State('custom-cnp-graph-purity', 'value'), 'Purity')
        self.add_autofill('cnp-plot-button', State('custom-cnp-graph-ploidy', 'value'), 'Ploidy')

    def set_default_review_data_annotations(self):
        self.add_review_data_annotation(
            annot_name='Purity',
            review_data_annotation=DataAnnotation('float', validate_input=validate_purity))
        self.add_review_data_annotation(
            annot_name='Ploidy',
            review_data_annotation=DataAnnotation('float', validate_input=validate_ploidy))

    def set_default_review_data_annotations_app_display(self):
        self.add_review_data_annotations_app_display(annot_name='Purity', app_display_type='number')
        self.add_review_data_annotations_app_display(annot_name='Ploidy', app_display_type='number')

        