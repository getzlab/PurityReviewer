"""
A reviewer dashboard that displays generic sample data and a allelic copy ratio profile 
for a given sample. The allelic copy ratio profile is linked to a table with the solutions 
from ABSOLUTE (Carter, 2014), where you can select a row and the corresponding ABSOLUTE "comb" 
solution will be plotted over the allelic copy ratio plot. You can filter the absolute solutions
displayed, by using a slider to select the range of solutions shown. If the absolute solution's 
purity value is within range of the a precalled purity value then the absolute solution is 
displayed in the absolute solution report. 
"""
from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from AnnoMate.ReviewerTemplate import ReviewerTemplate
from AnnoMate.AppComponents.DataTableComponents import gen_annotated_data_info_table_component
from AnnoMate.AnnotationDisplayComponent import NumberAnnotationDisplay, SelectAnnotationDisplay, TextAreaAnnotationDisplay

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

from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
import os
import pickle
from typing import Union, List, Dict
import sys

from PurityReviewer.AppComponents.AbsoluteSolutionsReportComponent import gen_absolute_solutions_report_component
from PurityReviewer.AppComponents.AbsoluteCustomSolutionComponent import gen_absolute_custom_solution_component
from PurityReviewer.AppComponents.AbsolutePrecalledSolutionComponent import gen_absolute_precalled_solutions_report_component

from PurityReviewer.AppComponents.utils import gen_cnp_figure, gen_mut_figure, parse_absolute_soln, validate_purity, validate_ploidy, CSIZE_DEFAULT


class PrecalledPurityReviewer(ReviewerTemplate):
    """
    Dashboard to allow user to view absolute solutions within a user-defined range of the 
    precalled purity value; default range is +/- 5% of the precalled purity value

    Carter, S., Cibulskis, K., Helman, E. et al. Absolute quantification of somatic DNA alterations in human cancer. Nat Biotechnol 30, 413–421 (2012). https://doi.org/10.1038/nbt.2203
    """
    def gen_data(self,
                 description: str,
                 df: pd.DataFrame,
                 index: List,
                 annot_df: pd.DataFrame = None,
                 annot_col_config_dict: Dict = None,
                 history_df: pd.DataFrame = None,
                 ) -> GenericData:
        """
        Generate data for PrecalledPurityReviewer object

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
                       sample_info_cols: List[str],
                       acs_col,
                       rdata_fn_col,
                       mut_fig_hover_data=None,
                       custom_parse_absolute_soln=None,
                       step_size=None
                       ) -> ReviewDataApp:
        """
        Parameters
        ==========
        sample_info_cols: list[str]
            List of columns in data
            
        acs_col: str
            Column name in data with path to seg file from alleliccapseg or other tsv with allelic copy ratio measurements

        rdata_fn_col: str
            Column name in data with LOCAL path to maf file. Should be predownloaded at set_review_data()
            
        mut_fig_hover_data: list
            List of column names to add to plotly hover_data in mutation figure
            
        custom_parse_absolute_soln: function
            Custom absolute parser function (rdata_path -> data_df)
            
        step_size: float
            Minimum increment allowed for purity (default is 0.01)

        Returns
        =======
        ReviewDataApp
            Review data app object to render the dashboard
        """
        app = ReviewDataApp()

        # add component that gets purity values within a specific range
        mut_fig_hover_data = [] if mut_fig_hover_data is None else mut_fig_hover_data
        app.add_component(
            gen_absolute_precalled_solutions_report_component(),
            rdata_fn_col=rdata_fn_col,
            acs_col=acs_col, 
            mut_fig_hover_data=mut_fig_hover_data,
            custom_parse_absolute_soln=custom_parse_absolute_soln,
        )

        # add component that displays a table with the annotated data info 
        app.add_component(
            gen_annotated_data_info_table_component(), 
            cols=sample_info_cols, 
            data_attribute='df',
            link_display_name=None
        )

        # add component that allows the user to create manual custom absolute solutions
        app.add_component(
            gen_absolute_custom_solution_component(step_size=step_size),
            acs_col=acs_col,
            step_size=step_size,
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
        self.add_annotation_display_component('Purity', NumberAnnotationDisplay())
        self.add_annotation_display_component('Ploidy', NumberAnnotationDisplay())
        self.add_annotation_display_component('Method', SelectAnnotationDisplay())
        self.add_annotation_display_component('Absolute_solution_idx', NumberAnnotationDisplay())
        self.add_annotation_display_component('Notes', TextAreaAnnotationDisplay())
