# secRecon: Secretory Pathway Reconstruction

This repository, **secRecon**, contains all the necessary code to generate the figures and analyses presented in the paper: *A reconstruction of the mammalian secretory pathway
identifies mechanisms regulating antibody production*

This paper introduces secRecon, a comprehensive reconstruction of the mammalian secretory pathway, highlighting its utility in contextualizing omics data and uncovering subprocess- to gene-level insights into secretory phenotypes. The repository includes multiple modules to facilitate the analysis of the secretory pathway over multiple species, cell types, and biological conditions. The repository is structured into several folders, each applying secRecon to a different dataset and addressing specific questions linking secretory machinery to phenotypes of interest.

### Supplementary Data
In addition to the code and scripts provided in this repository, the supplementary data accompanying the analyses are hosted on Synapse.org and can be accessed via the DOI: https://doi.org/10.7303/syn64026567.

## 1 - Network Reconstruction

This folder, **01 - Network Reconstruction**, contains all the essential code used to map data from different databases and datasets into secRecon and the Functional and PPI networkl topology generation.

### Contents

1. **1.1 - Feature\_Extraction.ipynb**: This Jupyter notebook performs ortholog mapping, extracts relevant features from different databases, and includes it into secRecon

2. **1.2 - Functional\_Ontology\_Network.ipynb**: This notebook is responsible for constructing the **Functional Ontology Network**. It uses the processed gene and process data to build a network that describes relationships between genes based on their biological function and ontology.&#x20;

3. **1.3 - PPI\_Ontology\_Network.ipynb**: The third notebook in the series focuses on integrating **Protein-Protein Interactions (PPI)** with ontology data to create a comprehensive **PPI Ontology Network**.&#x20;

## 2 - CHO vs Plasma

This folder contains notebooks and scripts used to perform various analyses on the CHO vs Plasma multi-omics dataset [1], aiming to explore differential secretory topologies and identify potential CHO engineering targets using secRecon.

### Contents

1. **2.1 - preprocess\_multiomics.ipynb**: This Jupyter notebook focuses on preprocessing the transcriptomics and proteomics datasets, grouping the data into upregulated and downregulated genes/proteins in human and murine plasma cells relative to antibody secreting CHO cells. This data is used dowmstream in network analysis (2.5)

2. **2.2 - correlation\_hclust.Rmd**: This R Markdown file performs Spearman correlation and dendrogram correlations of secretory pathway topologies at the gene/protein and process level for both transcriptomic and proteomic datasets. Gene Set Variation Analysis [2] (GSVA) is used to score secRecon geneset activity.

3. **2.3 - gsva\_limma.Rmd**: This R Markdown file performs differential analysis of secRecon GSVA scores using `limma` to identify differentially active and inactive secretory process in plasma cells vs CHO cells. Potential engineering targets in these enriched processes are identified.

5. **2.4 - Cho\_vs\_plasma\_network\_analysis.ipynb**: This Jupyter notebook generates and compares network features between CHO cells and plasma cells. It aims to highlight significant topological differences that could be leveraged to enhance protein production or understand cell-specific characteristics of the secretory pathway.

### Data
Multi-omics data used in this analysis was generated by Raab et al. and can be obtained from *10.1016/j.ymben.2024.03.007* Supplementary Materials; specific files required: *1-s2.0-S1096717624000521-mmc2.xlsx, 1-s2.0-S1096717624000521-mmc3.xlsx* 

[1] Raab, N. et al. Nature as blueprint: Global phenotype engineering of CHO production cells based on a multi-omics comparison with plasma cells. Metab. Eng. 83, 110–122 (2024). doi: 10.1016/j.ymben.2024.03.007

[2] Hänzelmann, S., Castelo, R. & Guinney, J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics 14, 1–15 (2013).

## 3 - SEC-seq Analysis

This folder, contains scripts used to perform analyses with a SEC-seq dataset [3] using secRecon to identify secretory pathway signatures linked to plasma cell differentiation and single-cell IgG secretion.

### Contents

1. **0a\_extract\_ppi\_network.Rmd**: This R Markdown file prepares the Protein-Protein Interaction (PPI) network edgelist for secRecon PPI scoring. PPI networks are extracted from STRINGdb and PCNet databases.

2. **0b\_Preprocess-DimensionReduction-Pseudotime.ipynb**: This Jupyter notebook preprocesses the raw SEC-seq data from 3 human donor samples, performs dimension reduction and pseudotime analysis as described by Cheng et al. 

3. **1\_score-secRecon-activity.ipynb**: This notebook scores the activity of secRecon processes using Seurat Feature scoring and ORIGINS2 PPI scoring [4]. secRecon ontology and annotations are used to defined genesets. Process scores are used as features for Dominance Analysis and Spearman correlation against single-cell phenotypes (IgG secretion, pseudotime).

4. **2\_Dominance-Analysis.ipynb**: This notebook performs Dominance Analysis to predict the most important secRecon processes linked to IgG secretion and plasma cell differentiation. Dominance Analysis [5] is applied to different clusters to dissect secretory pathway signatures associated with transcriptomic and phenotypic heterogeneity.

5. **3\_Correlation-DEG.ipynb**: This notebook performs secRecon subprocess- and gene-level correlation with IgG secretion and plasma cell differentiation phenotypes. Key secretory processes and machinery are identified as markers for IgG secretion efficiency and plasma cell differentiation.

6. **4\_IgGpopulation\_correlation.Rmd**: This R Markdown file performs secRecon gene-level correlation and dendrogram correlation for pseudobulk groups for IgG producing Leiden clusters. It aims to explore secretory topological differences between these distinct clusters.

### Data
SEC-seq data used in this analysis was generated by Cheng et al. and can be downloaded from GEO under accession number “GSE229042”. 

[3] Cheng, R. Y.-H. et al. SEC-seq: association of molecular signatures with antibody secretion in thousands of single human plasma cells. Nat. Commun. 14, 3567 (2023). doi: 10.1038/s41467-023-39367-8

[4] Senra, D., Guisoni, N. & Diambra, L. Cell annotation using scRNA-seq data: A protein-protein interaction network approach. MethodsX 10, 102179 (2023). (https://github.com/danielasenraoka/ORIGINS2)

[5] https://dominance-analysis.github.io/dominance-analysis/

## 4 - secRecon\_TF\_enrichment

This folder contains lists of transcription factors (TFs) enriched from IPA upstream regulator analysis, ChEA3, and Lund et al. RECON for secRecon genes. The aim is to identify transcription factors that may play a regulatory role in the mammalian secretory pathway and may be integrated/annotated in the knowledgebase in future versions.

### Contents

1. **secRecon\_TF\_enrichment.ipynb**: This Jupyter notebook compiles and compares these lists of TFs.

Krämer, A., Green, J., Pollard, J., Jr & Tugendreich, S. Causal analysis approaches in Ingenuity Pathway Analysis. Bioinformatics 30, 523–530 (2014).

Keenan, A. B. et al. ChEA3: transcription factor enrichment analysis by orthogonal omics integration. Nucleic Acids Res. 47, W212–W224 (2019).

Lund, A. M. et al. Network reconstruction of the mouse secretory pathway applied on CHO cell transcriptome data. BMC Syst. Biol. 11, 37 (2017).

## 5 - Gtex\_correlation

This folder contains scripts and GTEx data to explore potential systemic and literature biases in secRecon curation and relevance scoring.

### Contents

1. **0\_gtex\_secrecon\_char.Rmd**: This R Markdown file uses GTEx data to performs tissue-specific analysis, correlation average secRecon gene expression with curated metrics including number of annotated secRecon terms, max relevance score, mean relevance scoring.

### Data
GTEx RNA-seq data (TPM expression) can be downloaded from HPA database, consists of transcriptomics data spanning 35 tissues based on 46 tissue subtypes (https://www.proteinatlas.org/humanproteome/tissue/data#gtex_tissue_groups_rna)

GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318–1330 (2020).

GTEx Consortium et al. Genetic effects on gene expression across human tissues. Nature 550, 204–213 (2017).

## Usage Instructions

1. Clone this repository.
2. Navigate to the relevant folder for the analysis you wish to perform.
3. Make sure you have the appropriate software installed (e.g., Jupyter Notebook, R) along with the required Python and R packages. You can install the Python dependencies via:
   ```bash
   pip install -r requirements.txt
   ```
4. Run the Jupyter notebooks or R Markdown files in the specified order to ensure proper data processing and analysis.

## Unified Requirements

- Python 3.8+
- R version 4.0+
- Jupyter Notebook
- Required Python packages (install via `pip` or conda): `pandas`, `numpy`, `networkx`, `matplotlib`, `scikit-learn`, `scanpy`, `seaborn`, `scipy`, `pickle`, `dominance_analysis`
- Required R packages: `limma`, `GSVA`, `ggplot2`, `dplyr`, `corrplot`, `EnhancedVolcano`, `dendextend`

---
