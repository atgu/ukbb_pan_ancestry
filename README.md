# UKBB-Diverse-pops
GWAS across all analyzable phenotypes and diverse ancestry groups in the UK Biobank.

A full description of how these summary statistics were generated and how to work with them can be found here: http://pan.ukbb.broadinstitute.org/. 

# Code overview
The code in this repo can be used to generate all the figures in the manuscript. They look for data files in ./data, and if they are not found, downloads them as needed from the public data repository.

Each figure panel can be generated individually, or figures as a whole. For instance, in R/fig3_spectrum.R, we provide a function called `figure3()` which can generate the entirety of figure 3. Alternatively, running the code inside the function can generate each figure panel separately. Note that for some figures, on some R setups, attempting to generate the full figure by calling the function directly can crash R: we are uncertain of the cause of the issue, but it can be resolved by running the code inside the function step-wise.

# System Requirements
## Hardware requirements
There are no specific hardware requirements, except enough RAM to fit the datasets in memory (~16Gb).

## Software requirements
### OS Requirements
This package has been tested on macOS Sonoma (14.5).

### R Dependencies
The code depends on R packages, including those listed at the top of `constants.R`

# Installation Guide:

```
git clone https://github.com/atgu/ukbb_pan_ancestry
cd ukbb_pan_ancestry
R
```
- Install any packages needed, and run `constants.R` - if this runs without error, the remaining files will load.

## Citation

To cite the Pan-UKB resource, please use the citation: [Karczewski, Gupta, Kanai et al., 2024](https://medrxiv.org/content/10.1101/2024.03.13.24303864).
