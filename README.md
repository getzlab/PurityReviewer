# PurityReviewer

Suite of reviewers for reviewing purities.

![Purity Jupyter Reviewer](https://github.com/getzlab/JupyterReviewer/blob/master/images/ezgif.com-gif-maker.gif)

# Install

## Activate or Set up Conda Environment

This is **_highly_** recommended to manage different dependencies required by different reviewers.

See [Set up Conda Environment](https://github.com/getzlab/JupyterReviewer/blob/master/README.md#set-up-conda-environment) for details on how to download conda and configure an environment.
    
## Install PurityReviewer

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

# Custom and advanced usage

See `PurityReviewer/AppComponents` for pre-built components and their customizable parameters, and additional utility functions. 

For customizing annotations, adding new components, and other features, see [Intro_to_Reviewers.ipynb](https://github.com/getzlab/JupyterReviewer/blob/master/example_notebooks/Intro_to_Reviewers.ipynb).

For more detailed tutorials, see [AnnoMate](https://github.com/getzlab/AnnoMate)
