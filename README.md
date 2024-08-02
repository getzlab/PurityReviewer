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

## Why call purities?
The presence of normal cells in tumor samples dilutes the signal for tumor specific events such as mutations and copy number events. Estimating the purity (the percentage of cells in a tumor sample that are tumors) and ploidy (the average number of copies of the genome across tumor cells) allows us to adjust the raw sequencing data to infer relevant genomic characteristics of the tumor. These characteristics include (i) the absolute number of copies of each chromosome and (ii) clonal and subclonal architecture<sup>1</sup>. 

## Using ABSOLUTE<sup>1</sup> to call purities
ABSOLUTE<sup>1</sup> identifies likely purity and ploidy solutions for a tumor sample given information from bulk DNA sequencing data, including the allelic copy ratio (ACR) and variant allele fractions (VAF, proportion of reads supporting the mutation) for mutation calls. 

## Why manual review?
Different purity and subclonal architecture in a tumor sample can lead to similar observed mutation VAFs and copy number. For example, whole genome doubling events still correspond to “balanced” ACR, and look the same as a normal diploid genome. This ambiguity means there may be multiple purity-ploidy solutions with similar likelihood. ABSOLUTE may identify multiple likely purity-ploidy solutions, so an analyst has to decide among the solutions based on heuristics and prior biological/clinical knowledge about the tumor and sequencing artifacts. 

The *PurityReviewer* was created to automate this manual review process, allowing a reviewer to quickly render relevant complementary data and intermediate figures to help select the most appropriate purity solution. 

## Running ABSOLUTE
See [ABSOLUTE](https://github.com/getzlab/ABSOLUTE) for more details about input data. In brief:
1. Run a mutation calling tool or pipeline (e.g. Mutect1<sup>2</sup>). Perform additional artifact filtering so that the final mutation table contains only somatic mutations.
2. Run a copy number pipeline that produces a segmentation of the genome and calculates an ACR for each homolog at each segment (GATK CNV<sup>3</sup> or AllelicCapSeg<sup>4</sup>).
3. Run ABSOLUTE.

The [CGA characterization pipeline](https://github.com/broadinstitute/CGA_Production_Analysis_Pipeline) automatically calls mutations, produces the copy number profile, and runs ABSOLUTE. The output will be an RData table with the purity-ploidy solutions and corresponding annotations. 

## Reviewing purities

### General best practices when reviewing ABSOLUTE solutions 

#### Pick simplest solution
- Prefer solutions with integer CN that line up with "comb" peaks
- Largest peaks should be at reasonable CNs (e.g. 1 (haploid), 2 (diploid), 4 (tetraploid))
- Bottom peak should generally be close to zero (corresponding to clonal deletions)

#### Evidence from mutation multiplicities
- The multiplicity peak should line up with CN=1 
- There may possibly be a small peak at CN=2
- Purity ÷ 2 should align with the highest-VAF peak in balanced regions

#### Whole Genome Duplication (WGD) evidence
- WGD call should be supported by presence of clonal mutations (at least some at both CN = 1, 2 depending on timing of doubling)
- Use serial samples if you have them to resolve discrepancies (just remember varying possibilities for subclones)

#### Investigate/validate anomalies and discrepancies
- Check general QC, upstream tasks like AllelicCapSeg plots, etc.
- Look for orthogonal evidence in single cell data or serial samples
- (rare) A large subclonal homozygous deletion is almost certainly on sibling clones since cells cannot survive without at least 1 copy of large chromosomal regions
    - If sibling clones are unlikely (due to CCFs), investigate cause of artifact
- Both homologs subclonally deleted could be in separate clones

Using the *MatchedPurityReviewer*, you can automatically render each purity-ploidy solution at a time by clicking through an interactive table. The Reviewer includes an “autofill” button that you can click on to fill in your annotations directly once you have decided on a solution.

### Picking solutions manually 
If you have gone through all the ABSOLUTE solutions and have not found a purity that looks good for a certain sample, you can estimate your own purity-ploidy by deciding where the 0-line and 1-line of the comb plot should lie in the allelic copy ratio scale. Ideally, clonal deletion or LOH events correspond to the 0-line, and biallelic non-amplified, balanced segments correspond to the 1-line. 

To estimate the purity ($\alpha$):

$$\delta = l_1 - l_0$$

$$\alpha = 1 - \frac{l_0}{l_1}$$

To estimate average ploidy of the **whole sample** ($\tau$):

$$\tau = \frac{2 (1 - l_0)(1 - \alpha)}{\alpha l_0}$$

The *PurityReviewer* allows you to “draw” the 0-line and 1-line on top of an ACR plot, and automatically calculates the corresponding purity and ploidy. Once you have set the lines, the reviewer includes an autofill button to fill the annotation table with your custom solution.


### No solution
It is also possible that there is no clear/good solution for a given sample, especially if the sample’s ACR profile is too noisy (highly segmented) or the sample has very low purity. If this is the case, you should investigate the upstream data and tools for confirmation of the lack of signal; perhaps there was a bug in the pipeline that caused this issue. However, in the more likely event that the sample’s read count or purity is too low, you can decide to discard the sample.

The *PurityReviewer* includes a notes section to annotate this decision. Alternatively, you can add a custom annotation to label samples for specific action (e.g.: a checkbox with the option “no solution” or “discard”) and leave the remaining fields blank.

## PurityReviewer Annotations
The following is a list of the default annotations included in the *PurityReviewer*. Note that this reviewer has autofill buttons that will automatically fill in all of the following annotations (aside from notes) based on the solutions you have chosen in the reviewer as detailed in the previous section. 
- Purity: The percentage of tumor cells in a given sample 
- Ploidy: The average amount of DNA in the tumor cells (2 = diploid, 4 = tetraploid, etc.)
- Method: How you determined your solution - ABSOLUTE or Manually 
- ABSOLUTE solution index: The ABSOLUTE solution you chose (if you chose one)
- Notes: Any information to note about the sample

It is also possible to add your own annotations in addition to the default ones. The following is a list of additional annotations that could be useful to add depending on the specific project you are working on or the type of cancer you are studying:
- Whole genome doubling (yes/no)
- Quality (low/high)
- ChrX loss or gain

# Custom and advanced usage

See `PurityReviewer/AppComponents` for pre-built components and their customizable parameters, and additional utility functions. 

For customizing annotations, adding new components, and other features, see [Intro_to_AnnoMate_Reviewers.ipynb](https://github.com/getzlab/JupyterReviewer/blob/master/tutorial_notebooks/Intro_to_AnnoMate_Reviewers.ipynb).

For more detailed tutorials, see [AnnoMate](https://github.com/getzlab/AnnoMate)

# References
1. Carter, S. L. et al. Absolute quantification of somatic DNA alterations in human cancer. Nat. Biotechnol. 30, 413–21 (2012). doi: 10.1038/nbt.2203.
2. Cibulskis, K., Lawrence, M. S., Carter, S. L., Sivachenko, A., Jaffe, D., Sougnez, C., Gabriel, S., Meyerson, M., Lander, E. S., & Getz, G. (2013). Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nature biotechnology, 31(3), 213–219. https://doi.org/10.1038/nbt.2514
3. Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media
4. Carter, S., Meyerson, M. & Getz, G. Accurate estimation of homologue-specific DNA concentration-ratios in cancer samples allows long-range haplotyping. Nat Prec (2011). https://doi.org/10.1038/npre.2011.6494.1
