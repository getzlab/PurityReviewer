# PurityReviewer

Suite of reviewers for reviewing purities.

![Purity Jupyter Reviewer](https://github.com/getzlab/JupyterReviewer/blob/master/images/ezgif.com-gif-maker.gif)

# Install

## Activate or Set up Conda Environment

This is **_highly_** recommended to manage different dependencies required by different reviewers.

See [Set up Conda Environment](https://github.com/getzlab/JupyterReviewer/blob/master/README.md#set-up-conda-environment) for details on how to download conda and configure an environment.
    
## Install PurityReviewer

**Install R**

On ubuntu:
```
sudo apt update
sudo apt -y upgrade
sudo apt -y install r-base
```

Clone 
```
git clone git@github.com:getzlab/PurityReviewer.git 

# or in an existing repo
git submodule add git@github.com:getzlab/PurityReviewer.git
```

Install
```
cd PurityReviewer
conda activate <your_env>
pip install -e .
```

# Basic usage

See [purity_reviewer_example.ipynb](https://github.com/getzlab/PurityReviewer/tree/master/example_notebooks/purity_reviewer_example.ipynb) for basic examples and demos of the purity reviewers.

See `PurityReviewer/Reviewers` to see available pre-built reviewer options.

See `PurityReviewer/DataTypes` to see pre-built data configurations for purity review.

## Inputs

Sample dataframe `df`: index values correspond to samples or pairs, and must include a column with paths to allelic copy ratio segmentation tsv files. These files must have an AllelicCapSeg-like format. 

Each row corresponds to a segment.

- `Chromosome`: Chromosome of the segment
- `Start.bp`: Start position of the segment
- `End.bp`: End position of the segment
- `f`: median het site vaf
- `tau`: local total copy ratio
- `sigma.tau`: variance on tau
- `mu.major`: allelic copy ratio of the major allele
- `sigma.major`: variance on the major allele's allelic copy ratio
- `mu.minor`: allelic copy ratio of the minor allele
- `sigma.minor`: variance on the minor allele's allelic copy ratio

Possible segmentation tools include:
- AllelicCapSeg: https://github.com/aaronmck/CapSeg/blob/master/R/allelic/AllelicCapseg.R
- GATK ACNV: https://sites.google.com/a/broadinstitute.org/legacy-gatk-documentation/gatk-4-beta/7387-Description-and-examples-of-the-steps-in-the-ACNV-case-workflow

The Sample dataframe may also include other relevant data you may want to include during your review. Pass those columns as a list to `sample_info_cols` in `reviewer.set_review_app()` in your notebook. For each sample, the corresponding data in those columns will be displayed in a table in the dashboard.

# Custom and advanced usage

See `PurityReviewer/AppComponents` for pre-built components and their customizable parameters, and additional utility functions. 

For customizing annotations, adding new components, and other features, see [Intro_to_AnnoMate_Reviewers.ipynb](https://github.com/getzlab/JupyterReviewer/blob/master/tutorial_notebooks/Intro_to_AnnoMate_Reviewers.ipynb).

For more detailed tutorials, see [AnnoMate](https://github.com/getzlab/AnnoMate)
