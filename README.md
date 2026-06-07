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

1. **Covariate shift analysis.** For each gene, a linear model is fit to the retained loop patterns predicting that gene's logit-transformed loop weight from the sample-level statistics; coefficient t-values for each covariate are aggregated (e.g., all "age" t-values across loops) and the ordered gene statistics are input to GSEAs such as GO, KEGG or Reactome.

2. **Survival modelling.** Cox regression models of clinical covariates and gene weights identify genes significantly associated with clinical outcomes and prognosis.

In datasets where samples are dominated by a single loop (as was the case in our breast cancer study) our **Fourier enrichment** approach can detect which pathways/terms cluster in the dataset.

![Fourier enrichment workflow](figures/fourier.png)

- **Left:** A breast cancer dataset loop in a two-dimensional embedding, where each point represents a gene assigned an angular coordinate ($\theta_g$). Genes in the endometrial cancer pathway (highlighted in red) cluster near $\theta_g = 180^\circ$. 

- **Center-left and center-right:** Pathway membership and $\theta_g$ values are used to fit a Fourier model (intercept, $\cos\theta_g$, and $\sin\theta_g$ terms) and a baseline intercept-only model. A likelihood ratio test (LRT) compares the models, and resulting pathway p-values. 

- **Right:** Pathway p-values are corrected for multiple testing to identify Fourier-enriched pathways. For this loop, the endometrial cancer pathway is significantly enriched.

---

## Usage

The script `mola.R` contains code to run MOLA analyses.

### Requirements

- R >= 4.4.2
- Tested on Mac OSX and Linux systems - the maTilDA dependency doesn't seem to build on Windows.
- Significant RAM memory for larger numbers (> 500) of genes.

The top of `mola.R` installs (missing) necessary R packages and sets up the necessary reticulate conda environment - Python 3.13 must be installed on your system prior to running. Additionally, you need to clone the repository [maTilDA](https://github.com/IBM/matilda) into the same parent directory as MOLA.

### Clone Repository

```bash
git clone https://github.com/shaelebrown/MOLA.git
cd MOLA
```

---

## Quick Start

Here is a minimal example of how to identify loops, and their corresponding gene weight vectors, in several multiomic datasets. This code assumes that the directory '../multiomic_datasets' contains files sample1.csv, sample2.csv, ..., sample50.csv, where all rows correspond to the same set (and order) of genes and the columns correspond to the same omic measures (e.g., "RNAseq", "ATACseq", etc.).

```r
source("mola.R")

# process all samples
lapply(X = 1:50, FUN = function(X){

  # read in that sample's data
  # the genes can be the rownames, but should not be in a column "X"
  df <- read.csv(paste0("../multiomic_datasets/sample", X, ".csv"))
  
  # assume directory "../mola_results" already exists
  results_dir <- paste0("../mola_results/sample", X)
  
  # run MOLA and save results
  # keep all loops to jointly filter across all samples
  run_mola(df = df, results_dir = results_dir, alpha = 1.0, verbose = T)

})

```

Results for sample `X` will be saved in the directory "../mola_results/sample`X`", with the following files:

- **gene_weights.RData** - the gene weights for each identified (and retained if $alpha < 1$) loop. Loading this file will create an object called `point_weights_mat`, in which rows are the genes (in the same row order as all the samples), columns are the loop (named by loop key) and entries are the weights. 
- **gene_interaction_weights.RData** - the gene interaction (i.e., edge) weights for each loop.
- **persistence_diagram.csv** - the computed persistence diagram with loop key as a column.
- **params.RData** - the alpha and enclosing radius used in calculations.
- **log.txt** - informative messages and diagnostics in the case of errors.

To filter loops jointly across all samples, you can use the `filter_loops` function, like so:

```r
# read in the diagrams for all samples
joint_diagram <- lapply(X = 1:50, FUN = function(X){

  diag <- read.csv(paste0('../mola_results/sample', X, '/persistence_diagram.csv'))
  diag$samp <- X
  return(diag)

})

# filter loops - subsets joint_diagram for loop p-values < 0.05
filtered_diagram <- filter_loops(diagram = joint_diagram, alpha = 0.05)

# the columns "samp" and "key" will give the necessary identifiers for
# using retained loops
weights <- do.call(rbind, lapply(X = 1:nrow(filtered_diagram), FUN = function(X){

  samp <- filtered_diagram$samp[[X]]
  key <- filtered_diagram$key[[X]]
  
  # read in that sample's point_weights_mat and subset for the key
  load(paste0("../mola_results/sample", samp, "/gene_weights.RData"))
  harmonic_weights <- point_weights_mat[, as.character(key)]
  harmonic_weights <- harmonic_weights/max(harmonic_weights)
  return(matrix(data = harmonic_weights, columns = 1))

}))

```

For covariate shift analysis you will need to modify the *covariate_shift_analysis* function to adjust for your own covariates and modelling procedures - outputs will be a file **tvals.RData** stored in a user-supplied results directory:

```r
# read in sample-level statistics file
pheno <- read.csv('../pheno.csv')

# match with retained loops, assuming "samp" is the sample ID column in
# both dataframes
pheno <- merge(filtered_loops, pheno, by = c("samp"))

# read in or define gene vector, corresponding to the rows of "weights"
genes <- c("ERLIN2", "CLIC6", "SCGB2A2", "GSTM1", "ZNF703", 
           "CDC42BPA", "ANKRD30A", "HSPA12B", "ASH2L", "IKZF1", 
           "VSTM2A", "BAG4", "PSMG3", "TACC1", "TAT", 
           "EPO", "FGF3", "TFAP2B", "SHH", "CCND1", 
           "ADIPOQ", "PCDHGB6", "ANO1", "CEACAM5", "MUC5B", 
           "HOXB2", "MT1G", "CALML5", "STARD3", "PGAP3", 
           "TCAP", "HOXB13", "HSPB6", "GRIA2", "IKZF3", 
           "NRG1", "GRB7", "CRABP1", "BMPR1B", "WSCD2", 
           "MYEOV", "BEST4", "GFRA1", "CDK12", "PSMD3", 
           "PCBP1", "A2ML1", "CLCA2", "GIMAP7", "FBXL20")
           
# define results directory to save the output file (results_dir/tvals.RData)
results_dir <- "../mola_results"
           
# run the covariate shift analysis
covariate_shift_analysis(weights = weights, pheno = pheno, genes = genes)

# load the results
load(paste0(results_dir, '/tvals.RData'))

# variables are named lists of gene t-values, each in decreasing order
age_tvals[1:5] # these variables depend on your modifications

```

For a survival analysis, you will need to modify the *survival_analysis* function to adjust for your own covariates and modelling procedures - outputs will be a list, either empty if no genes are significant after FDR correction (at the $\alpha = 0.05$ level) or named for significant genes as the Cox model results for those genes.

```r
# assuming pheno contains survival data columns PFI and PFI.time:
# (but adjust for your own variables)
significant_genes <- survival_analysis(surv_data = pheno, weights = weights, genes = genes)

# inspect
names(significant_genes)
significant_genes[["RET"]]

```

---

## Data Availability

The MOLA manuscript provided two case study analyses of an Influenza A virus (IAV) dataset and a TCGA breast cancer dataset.

The raw data from the IAV dataset are available from the European Genome-phenome Archive (EGA), under accession numbers EGAD00001008422 and EGAD00001008359.

The TCGA breast cancer dataset analyzed during the current study were accessed via the R package UCSCXenaTools version 1.7.0 (Wang, 2019) on March 12, 2026.

---

## Citation

If you use MOLA, please cite:

```bibtex
@article{MOLA2026,
 title   = {MOLA: a Novel Topological Data Analysis Framework for Analyzing Multiomic Loops in Precision Medicine},
  author  = {Brown, Shael and Xiao, Bowei and Oros Klein, Kathleen and Tonin, Patricia N. and Dupuis, Josée and Zhang, Qihuang and Greenwood, Celia M.T.},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.1101/2026.xx.xx.xxxxxx},
  note    = {Preprint}
}
```

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Contact

For questions please contact:

- Shael Brown
- shaelebrown@gmail.com

---

## Acknowledgments

The authors acknowledge funding from the DNA-to-RNA (D2R) Foundational Projects (Cycle 1) grant: Modelling individual differences in patterns of epigenetic regulation.

---

## Version History

### v1.0.0

- Initial release accompanying manuscript publication.

---

## References

Gurnari, D., Guzmán-Sáenz, A., Utro, F. et al. Probing omics data via harmonic persistent homology. Sci Rep 15, 38836 (2025).

Bobrowski, O., Skraba, P. A universal null-distribution for topological data analysis. Sci Rep 13, 12274 (2023).

Wang et al., (2019). The UCSCXenaTools R package: a toolkit for accessing genomics data from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq. Journal of Open Source Software, 4(40), 1627, https://doi.org/10.21105/joss.01627