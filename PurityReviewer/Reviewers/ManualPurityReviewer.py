"""
A reviewer dashboard that displays generic sample data and a allelic copy ratio profile for a given sample. The allelic copy ratio profile includes sliders and inputs to manually set the 0 and 1 line corresponding to the integer assignment of the genomic segments and automatically calculates a purity. 
"""
from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from AnnoMate.AppComponents.DataTableComponents import gen_annotated_data_info_table_component
from AnnoMate.ReviewerTemplate import ReviewerTemplate
from AnnoMate.AnnotationDisplayComponent import NumberAnnotationDisplay, TextAreaAnnotationDisplay

import pandas as pd
import numpy as np

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

from cnv_suite.visualize import plot_acr_interactive

from rpy2.robjects import r, pandas2ri
import os
import pickle

from typing import List, Dict

from dash.dependencies import State

class ManualPurityReviewer(ReviewerTemplate):
    """
    Dashboard to iterate through a sample/pair and manually select and annotate the best fit purity and ploidy solution. Either set the 0 and 1 line manually, or input a purity and ploidy.
    """
    def gen_data(self,
                 description: str,
                 df: pd.DataFrame,
                 index: List,
                 annot_df: pd.DataFrame = None,
                 annot_col_config_dict: Dict = None,
                 history_df: pd.DataFrame = None) -> GenericData:
        """Generate data for PurityReviewer object
        Parameters
        ==========
        df: pd.DataFrame
            Pandas dataframe with samples/pairs for rows, and relevant data in columns. 
            Columns should include an allelic copy ratio segmentation file (ie GATK ACNV pipeline alleliccapseg)

        Returns
        =======
        GenericData
            Data object that contains only one dataframe
        """
        pandas2ri.activate()

        return GenericData(index=index,
                           description=description,
                           df=df,
                           annot_df=annot_df,
                           annot_col_config_dict=annot_col_config_dict,
                           history_df=history_df)

    def gen_review_app(self,
                       sample_info_cols,
                       acs_col,
                       csize=None,
                       step_size=None
                       ) -> ReviewDataApp:
        """
        Parameters
        ==========
        sample_info_cols: list[str]
            List of columns in to display from data
            
        acs_col: str
            Column name in data with path to seg file from alleliccapseg or other tsv with allelic copy ratio measurements
            
        csize: dict
            Dictionary with chromosome sizes (see AppComponents.utils.CSIZE_DEFAULT for hg19)
            
        step_size: float
            Minimum increment allowed for purity (default is 0.01)

        Returns
        =======
        ReviewDataApp
            Review data app object to render the dashboard
        """
                           
        app = ReviewDataApp()
        
        app.add_component(
            gen_annotated_data_info_table_component(), 
            cols=sample_info_cols, 
            data_attribute='df',
            link_display_name=None
        )

        app.add_component(
            gen_absolute_custom_solution_component(step_size=step_size),
            acs_col=acs_col,
            step_size=step_size,
            csize=csize
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
        self.add_review_data_annotation(
            annot_name='Notes', 
            review_data_annotation=DataAnnotation(
                annot_value_type='string')
        )

    def set_default_review_data_annotations_app_display(self):
        self.add_annotation_display_component('Purity', NumberAnnotationDisplay())
        self.add_annotation_display_component('Ploidy', NumberAnnotationDisplay())
        self.add_annotation_display_component('Notes', TextAreaAnnotationDisplay())
        