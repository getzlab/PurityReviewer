import pandas as pd
from typing import List, Dict

from dash.dependencies import State
from rpy2.robjects import pandas2ri

from PurityReviewers.AppComponents.AbsoluteCustomSolutionComponent import gen_absolute_custom_solution_component
from PurityReviewers.AppComponents.utils import validate_purity, validate_ploidy
from JupyterReviewer.Data import DataAnnotation
from JupyterReviewer.ReviewDataApp import ReviewDataApp
from JupyterReviewer.DataTypes.GenericData import GenericData
from JupyterReviewer.ReviewerTemplate import ReviewerTemplate
from JupyterReviewer.AppComponents.DataTableComponents import gen_annotated_data_info_table_component
from JupyterReviewer.AnnotationDisplayComponent import NumberAnnotationDisplay


class ManualPurityReviewer(ReviewerTemplate):
    def gen_data(self,
                 description: str,
                 df: pd.DataFrame,
                 index: List,
                 annot_df: pd.DataFrame = None,
                 annot_col_config_dict: Dict = None,
                 history_df: pd.DataFrame = None) -> GenericData:
        """Generate data for PurityReviewer object

        Returns
        -------
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
            List of columns in data
        acs_col: str
            Column name in data with path to seg file from alleliccapseg
        csize: dict
            Dictionary with chromosome sizes
        step_size: float
            Minimum increment allowed for purity (default is 0.01)
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

    def set_default_review_data_annotations_app_display(self):
        self.add_annotation_display_component('Purity', NumberAnnotationDisplay())
        self.add_annotation_display_component('Ploidy', NumberAnnotationDisplay())
        