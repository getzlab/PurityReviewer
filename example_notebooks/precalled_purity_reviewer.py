import os
import sys

sys.path.append(os.path.abspath( os.path.join("..") ))

from PurityReviewer.Reviewers.PrecalledPurityReviewer import PrecalledPurityReviewer
from PurityReviewer.AppComponents.utils import CSIZE_DEFAULT
from PurityReviewer.AppComponents.utils import download_rdata
from PurityReviewer.AppComponents.utils import add_precalled_purities_to_pairs

import dalmatian
import pandas as pd
import numpy as np
import os

# Uncomment the following line to set the google project if you haven't already done so 
# os.environ["GCLOUD_PROJECT"] = <google project>

# get data from terra
workspace_name = 'TCGA_BRCA_WGS' # set this to the correct workspace name
workspace = f'broad-tcga-wgs-terra/{workspace_name}'
wm = dalmatian.WorkspaceManager(workspace)
wm_pairs_df = wm.get_pairs()
wm_samples_df = wm.get_samples()
wm_pairs_df = wm_pairs_df[wm_pairs_df['absolute_rdata_WGS'].notna()]

# download rdata locally
local_rdata_dir = 'local_rdata'
wm_pairs_df['local_absolute_rdata'] = download_rdata(wm_pairs_df['absolute_rdata_WGS'], rdata_dir=local_rdata_dir)

# set custom csize for sex chromosomes
sex_chr_map = {'23': 'X', '24': 'Y'}
rename_chroms = {x: sex_chr_map[x] if x in sex_chr_map.keys() else x for x in CSIZE_DEFAULT.keys()}
custom_csize = {f'chr{rename_chroms[chrom]}': length for chrom, length in CSIZE_DEFAULT.items()}
custom_csize

# associate precalled purities with pairs
wm_pairs_df = add_precalled_purities_to_pairs(wm_pairs_df, wm_samples_df)
wm_pairs_df = wm_pairs_df.set_index('tumor_submitter_id')

# create a reviewer instance and set the review data
precalled_purity_reviewer = PrecalledPurityReviewer()
precalled_purity_reviewer.set_review_data(
    data_path = 'precalled_purity_reviewer_output', 
    description= 'BRCA purity review', 
    df=wm_pairs_df,
    index=wm_pairs_df.index,
)

# specify column names if they are not the default
precalled_purity_reviewer.set_review_app(
    sample_info_cols=['absolute_highres_plot_WGS', 'hapaseg_allelic_segmentation_plot_WGS'],
    acs_col='hapaseg_segfile_WGS', 
    rdata_fn_col='local_absolute_rdata',
    mut_fig_hover_data=['Hugo_Symbol', 'Chromosome', 'Start_position'],
)

# set default annotations and autofill
precalled_purity_reviewer.set_default_review_data_annotations_configuration()
precalled_purity_reviewer.set_default_autofill()

# run the review app
precalled_purity_reviewer.run(port=8098)