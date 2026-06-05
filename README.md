# MOLA
### A novel topological data analysis framework for analyzing multiomic loops in precision medicine

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Repository accompanying the manuscript:

> Brown S, Xiao B, Klein K, Tonin P, Dupuis J, Zhang Q and Greenwood G. "MOLA: a novel topological data analysis framework for analyzing multiomic loops in precision medicine." bioRxiv, 2026.

---

## Overview

### What is MOLA?

A previous multiomic analysis framework (Gurnari, 2025) demonstrated the utility of *topological data analysis* (TDA) for identifying subtle non-linear structure in high-dimensional omic data, i.e., coordinated loop patterns of gene omic profiles represented as gene-valued vectors. However, methods for interpreting and relating these loop features to sample-level covariates remain limited.

To address this gap, we developed **MOLA** (**M**ulti**O**mic **L**oop **A**nalysis), a framework for the characterization and interpretation of multiomic loops. MOLA extends previous workflows with tools for:

* **Loop visualization** for exploratory analysis and interpretation
* **Loop filtering** to prioritize well-sampled structures
* **Covariate shift analysis** to quantify associations between loops and sample-level variables
* **Loop Fourier enrichment analysis** to identify pathways that cluster around loops

MOLA is designed for the analysis of richly-annotated multiomic datasets that include sample-level phenotypic, clinical, or experimental metadata. By linking topological structure to biological and clinical variables, MOLA enables the discovery of complex signals that may be overlooked by conventional linear approaches, with anticipated applications in precision medicine, biomarker discovery, and risk modeling.

---

## MOLA Workflow

![MOLA workflow](figures/mola_overview.png)

The initial steps of all MOLA analyses are:

(A) **Identification and visualization of loop structures.** *Left:* a single multiomic sample, projected onto three (RNA-seq, promoter methylation, and enhancer methylation) of the total 13 omic dimensions in the dataset, with each point representing a gene. Using the full 13-dimensional dataset, a gene–gene correlation distance matrix is computed (*center-left*) and provided as input to the Gurnari framework, which generates persistence diagrams and associated (harmonic) loop weight vectors. In the persistence diagram (*center*), each point represents a loop: the x-axis denotes the birth radius, the correlation distance at which genes connect around the loop circumference, while the y-axis denotes the death radius, the distance at which genes connect across the loop and fill it in. For a highlighted loop, the multiomic sample can be coloured by the Gurnari loop gene vector (*center-right*), and genes can be projected onto a 2D embedding of the loop (*right*). 

(B) **Loop filtering procedure.** A transformed loop statistic based on the death-to-birth ratio is computed for all loops across all samples (*left*), the distribution of which follows a LGumbel distribution (*right*, see Bobrowski, 2023). Only loops with unusually large death/birth ratios, including the highlighted example, are retained for downstream analysis.

Subsequent steps are flexible and depend on the study goals and available data; the MOLA paper includes two case studies showing examples of:

1. **Covariate shift analysis.** For each gene, a linear model is fit predicting that gene's logit-transformed loop weight across retained patterns from the sample-level statistics of the sample of each retained loop; coefficient t-values for each covariate are aggregated (e.g., all "age" t-values across loops) and the ordered gene statistics are input to GSEAs such as GO, KEGG or Reactome.

2. **Survival modelling.** Cox regression models of clinical covariates and gene weights identify genes associated with clinical outcomes and prognosis.

In datasets where samples are dominated by a single loop (as was the case in our breast cancer study) our **Fourier enrichment** approach can detect which pathways/terms cluster in the dataset.

---

## Usage

The script *mola.R* contains code to run MOLA analyses.

### Requirements

- R >= 4.4.2
- Tested on Mac OSX and Linux systems - maTilDA dependency doesn't seem to build on Windows.
- Significant RAM memory for larger numbers of genes.

The top of *mola.R* contains code to install necessary R packages and set up the necessary reticulate conda environment - Python 3.13 must be installed on your system prior to running. Additionally, you need to clone the repository [maTilDA](https://github.com/IBM/matilda) into the same parent directory as this directory.

### Clone Repository

```bash
git clone https://github.com/shaelebrown/MOLA.git
cd MOLA
```

---

## Quick Start

Here is a minimal example of how to identify loops, and their corresponding gene weight vectors, in a multiomic dataset:

```r
source("mola.R")

df <- read.csv('gene_omic_dataset.csv')
results_dir = '../mola_results/sample1'
alpha <- 0.05

run_mola(df = df, results_dir = results_dir, alpha = alpha, verbose = T)
```

Results will be saved in the directory "../mola_results/sample1", with the following files:

- **gene_weights.RData** - the gene weights for each identified (and retained) loop.
- **gene_interaction_weights.RData** - the gene interaction (i.e., edge) weights for each loop.
- **persistence_diagram.csv** - the computed persistence diagram.
- **params.RData** - the alpha and enclosing radius used in calculations.
- **log.txt** - informative messages and diagnostics in the case of errors.

For covariate shift analysis you will need to modify the *covariate_shift_analysis* function to adjust for your own covariates and modelling procedures - outputs will be a file **tvals.RData** stored in a user-supplied results directory:

```r
# FILL IN
run_mola(df = df, results_dir = results_dir, alpha = alpha, verbose = T)
```

For a survival analysis, WHAT

---

## Input Data

### Data Availability

Explain:

- Public datasets used
- Accessions
- Download links
- Controlled-access data restrictions

Example:

| Dataset | Source | Accession |
|----------|----------|----------|
| TCGA | NCI | XXXXX |
| GEO | NCBI | GSEXXXXX |

### Data Format

Describe expected input format.

Example:

```csv
sample_id,gene,expression,methylation,...
```

---

## Output Files

Describe all outputs generated by MOLA.

| Output | Description |
|----------|-------------|
| loops.csv | Detected loops |
| persistence.csv | Persistence statistics |
| network.rds | Loop network object |
| figures/ | Generated plots |

---

## Citation

If you use MOLA, please cite:

```bibtex
@article{MOLA202X,
 ...
}
```

DOI:

```
doi:xxxx
```

---

## License

Specify license.

Example:

MIT License

or

GNU GPL v3

---

## Contact

For questions:

- Name
- Institution
- Email

---

## Acknowledgments

Funding sources:

- NIH
- NSF
- Institutional support

Collaborators and contributors.

---

## Version History

### v1.0.0

- Initial release accompanying manuscript publication.

---

## References

Gurnari, D., Guzmán-Sáenz, A., Utro, F. et al. Probing omics data via harmonic persistent homology. Sci Rep 15, 38836 (2025).

Bobrowski, O., Skraba, P. A universal null-distribution for topological data analysis. Sci Rep 13, 12274 (2023).