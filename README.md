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

## Output Annotations

Note that this reviewer has autofill buttons that will automatically fill in all of the following annotations (aside from notes) based on the solutions you have chosen in the reviewer:
- Purity: The percentage of tumor cells in a given sample (Normal “contamination” produces shift in ACR)
- Ploidy: The average amount of DNA in the tumor cells (2 = diploid, 4 = tetraploid, etc.)
- Method: How you determined your solution - Absolute or manual 
- ABSOLUTE solution: The ABSOLUTE solution you chose (if you chose one)
- Notes: Any information to note about the sample 

It is also possible to add your own annotations in addition to the default ones. The following is a list of additional annotations that could be useful to add depending on the specific project you are working on / type of cancer you are studying:
- Whole genome doubling (yes/no)
- Quality (low/high)
- ChrX loss or gain

# Review With PurityReviewer

## Why is this type of review important 
Nearly all bulk tumor samples, tissue or blood, will contain some amount of normal cells. Subsequently, DNA or RNA collected from these tissues will be a mixture of both tumor and normal DNA/RNA sequences. Estimating the percentage of DNA that actually comes from the tumor (purity) is critical for downstream analysis and interpretation. For example, a mutation with a low (<0.2) variant allele fraction (VAF, percent of reads supporting the alternate allele) could really be in all the tumor cells if the sample also has low purity. Similarly, quantifying a purity will help determine whether copy number events (amplifications and deletions) are clonal (in all tumor cells) or subclonal (in only a subset). 

One of the tools in the CGA pipeline, ABSOLUTE (Carter et. al), was created to infer purities, along with other metrics, of bulk tumor samples that best explains the input data. This input data consists of somatic mutations (SNVs) and their VAFs, and segmentations of the genome using changes in total copy ratio (TCR) and allelic copy ratio (ACR) that signify copy number events. Since these ratios are relative, there may be multiple purity solutions that could produce the observed SNVs and Copy number events. Thus, ABSOLUTE returns multiple possible purity solutions, and manual review is required to select the most appropriate solution based on prior biological and/or clinical knowledge about the tumor and sequencing artifacts. The PurityReviewer was created to automate this manual review process, allowing a reviewer to quickly render relevant complementary data and intermediate figures to help select the most appropriate purity solution. 

## How is this review conducted 

### Pre-reviewer use: Pipelines and tools run before purity can be reviewed 
The steps preceding this review are: 
1. QC sequencing results: looking over initial BAM QC metrics to filter out any samples with low QC scores from downstream analyses
2. Check for sample swaps: fingerprinting samples within and across patients to identify potential mismatches during sequencing.
3. Run the CGA pipeline: the pipeline will automatically run ABSOLUTE and generate possible purity solutions to be reviewed, in addition to mutation calling and filtering. The output will be an RData table with the solutions and corresponding annotations. 

### General best practices reviewing ABSOLUTE solutions 

#### Pick simplest solution
1. Prefer solutions with integer CN that line up with comb peaks
2. Largest peaks should be at “easy” CNs (e.g. 1, 2, 4)
3. Bottom peak should generally be close to zero

#### Multiplicity peack should line up with CN=1 (possibly small peak at CN=2)
1. This plot is biased due to clonal/subclonal boundary - better to plot yourself if unclear
2. Similarly, ɑ/2 should divide highest-AF peak

#### Whole Genome Doubling/Trippling (WGD/WGT, etc.) should have support from clonal mutations (at least some at both CN = 1, 2 depending on timing of doubling)
1. Use serial samples if you have them to resolve discrepancies (just remember varying possibilities for subclones)

#### Investigate/validate anomalies and discrepancies
1. Check QC, upstream tasks like AllelicCapSeg plots, etc.
2. (rare) Large subclonal homozygous deletion - almost certainly on sibling clones 
    a. Can’t survive without at least 1 copy

### Using the PurityReviewer: How to use the reviewer, what to look for, tips and tricks 
Once you have all the information you need to review sample purities, you can run the PurityReviewer dashboard following the instructions in the README. This dashboard is where you can then render the data for each sample, select a purity solution, and make any additional annotations or notes. 

#### Picking ABSOLUTE solutions
Ideally, one of the ABSOLUTE solutions will be good and that will be what you pick for the purity. Since this is the goal, the first step of the purity review process is to go through each ABSOLUTE solution one by one, look at the corresponding copy number profiles and mutations, and determine if the solution looks good or not. In general, you want to select the solution where the “comb plot” - or the integer CN lines - match up with the most prominent segment cluster peaks and generally select the simplest solution that explains the data. If you’re not sure about a certain solution, you can look further at the linked AllelicCapSeg plots to see if the total copy ratios make sense with the allelic copy number in the main copy number plot on the dashboard. If you find a solution that looks good, you can simply click the “Pick current solution” button in the annotations panel to fill in all of the information for that solution, submit your annotation, and move on to the next sample. 

#### Picking solutions manually 
If after you have gone through all the ABSOLUTE solutions and have not found a purity that looks good, you can estimate your own purity in the “manual purity” section of PurityReviewer. To do so, you manually set the 0 line and 1 line to correspond to where you expect a clonal deletion and normal segments should lie in the allelic copy ratio scale, respectively. Once you estimate the purity this way, click the “Use manual solution” button in the annotations panel to fill in the information from your manual solution and submit your annotation. 

#### No solution
It is also possible that there is no clear/good solution for a given sample, especially if the sample has very low purity. If this is the case, you would likely want to discard the sample. You could add your own annotations (e.g.: a checkbox with the option “no solution” or “discard”), or simply write ‘discard’ in the “notes” section of the annotations panel. If this is the case, all other annotation fields can be left blank. 

## Next steps after review: Now that review has been conducted, what is done with the information found
Once you have gone through all of your samples, you can find all of your annotations in the annotations_df. This pandas DataFrame can then be merged with all of the original data in the Terra pairs table and uploaded back to the pairs table. This can either be done manually or with the use of dalmatian, a python package that allows for interaction with Terra directly through pandas DataFrames. 

Now that the original data is annotated with finding from your purity review, downstream processes that require purity and ploidy, like pre-phylogic and PhylogicNDT, can be run, filtering out the samples where low or unclear purity was found. This filtering can easily be done by looking for the ‘discard’ or ‘remove’ annotation. 

# Custom and advanced usage

See `PurityReviewer/AppComponents` for pre-built components and their customizable parameters, and additional utility functions. 

For customizing annotations, adding new components, and other features, see [Intro_to_AnnoMate_Reviewers.ipynb](https://github.com/getzlab/JupyterReviewer/blob/master/tutorial_notebooks/Intro_to_AnnoMate_Reviewers.ipynb).

For more detailed tutorials, see [AnnoMate](https://github.com/getzlab/AnnoMate)

# References
ABSOLUTE: Carter, S. L. et al. Absolute quantification of somatic DNA alterations in human cancer. Nat. Biotechnol. 30, 413–21 (2012). doi: 10.1038/nbt.2203.