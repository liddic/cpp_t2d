# cpp_t2d
Examining shared trends in compound-associated functional capacities of degraded ecosystem soil microbiomes and type 2 diabetes gut microbiomes

## About this repository

**Contains code for the article:**

*Shared potential metabolism trends in degraded soils and type 2 diabetes gut microbiomes*

by Craig Liddicoat, Bart A. Eijkelkamp, Timothy R. Cavagnaro, Jake M. Robinson, Kiri Joy Wallace, Andrew D. Barnes, Garth Harmsworth, Damien J. Keating, Robert A. Edwards and Martin F. Breed

Preprint server [to be updated]: [BIORXIV/2025/642605](https://doi.org/10.1101/2025.03.11.642605)

This repository documents code used in the above study, and is not intended for ongoing code development. Sub-folders contain scripts for bioinformatic processing and analysis of the case study datasets described below.

## Data sources
Datasets comprise soil metagenomes from five ecosystem quality gradients, and gut metagenomes from two cohorts of untreated T2D versus normal healthy subjects

**Soil datasets**
- Post-mining ecosystem restoration (data published by [Sun & Badgley 2019](https://www.sciencedirect.com/science/article/abs/pii/S0038071719301385?via%3Dihub))
- Disturbed versus natural soils (data from the [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/australian-microbiome) previously examined by [Liddicoat et al 2024](https://doi.org/10.1016/j.scitotenv.2024.173543)
- Prairie restoration (data published by [Mason et al 2023](https://academic.oup.com/femsec/article/99/11/fiad120/7288652))
- Plantation succession (data published by [Zuo et al 2026](https://link.springer.com/article/10.1007/s11104-025-08246-0))
- Vegetation succession (data published by [Wang et al 2024](https://academic.oup.com/ismecommun/article/4/1/ycae128/7833432))

**T2D datasets**
- Chinese cohort (data published by [Forslund et al 2015](https://www.nature.com/articles/nature15766) and [Qin et al 2012](https://www.nature.com/articles/nature11450))
- Swedish cohort (data published by [Forslund et al 2015](https://www.nature.com/articles/nature15766) and [Karlsson et al 2013](https://www.nature.com/articles/nature12198))

Study datasets are available from NCBI [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/) accessions PRJEB1786, PRJNA422434, PRJNA1215775, PRJNA1215778, PRJNA1215780, PRJNA1215781, PRJNA1080685; [MG-RAST](https://www.mg-rast.org/) project mgp16379; [JGI IMG](https://img.jgi.doe.gov/cgi-bin/m/main.cgi) Study ID Gs0144357, and the [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/australian-microbiome).

**Supplementary datasets**

Additional datasets were analysed to provide further validation of the new compound processing potential (CPP) method:
- Mice gut metagenomes raised in soil from desert, grassland and forest (data published by [Liu et al 2021](https://doi.org/10.1080/19490976.2020.1830699); NCBI Sequence read archive, accession: PRJNA542998)
- Soil metagenomes at 0h, 8h, 24h, and 48h following glucose amendment (data published by [Chuckran et al 2020](https://journals.asm.org/doi/10.1128/mra.00895-20); NCBI Sequence Read Archive accessions PRJNA539715, PRJNA539712, PRJNA539720, PRJNA539713, PRJNA539717, PRJNA539718, PRJNA539719, PRJNA539721, PRJNA539722, PRJNA539723, PRJNA539714, PRJNA539711, PRJNA539716.)
- Single-species gut bacteria culture metagenome standards (data published by [Amos et al 2020](https://link.springer.com/article/10.1186/s40168-020-00856-3) for the National Institute for Biological Standards and Control; NCBI Sequence Read Archive accession PRJNA622674)

## Method overview

**Compound processing potential (CPP)**

This study uses microbial metagenome functional profiling at the resolution of individual compounds. 

The CPP method (Fig. 1) uses outputs from [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS) functional profiling, which estimates the relative distribution of functional capacity across various biochemical pathways. Using associations defined in the [ModelSEED biochemistry database](https://github.com/ModelSEED/ModelSEEDDatabase/tree/master/Biochemistry), functional pathway relative abundances are divided among all linked reactions and compounds with weighting to account for molar ratios (stoichiometry), without thermodynamic or directional constraints.  Reaction inputs and products are treated equally, because often microbes can facilitate a string of reactions (so products become inputs, and so on). For all carbon (C)-containing compounds, elemental ratios of oxygen (O):C, hydrogen (H):C, and nitrogen (N):C are calculated to enable optional '3-d' visualisation and analysis of energetically and chemically similar compounds (Fig. 2). 

CPP values represent the hypothetical functional capacity (%) allocated to each compound, reflecting their potential metabolism by a given metagenome.

![cpp3d method overview](/ancillary_files/Figure1-CPP-workflow-Jun2026.png)

_**Figure 1.** Method overview_

![Example analyses](/ancillary_files/Figure2-Example-analyses.png)

_**Figure 2.** Example CPP analyses from soil microbiomes under post-mining forest ecosystem restoration (raw metagenome data from Sun & Badgley 2019). a) Potential metabolism for CO2 increases, suggesting rising microbial activity with revegetation age. b) Visualisation of potential metabolism (log10 CPP values) for 7,736 carbon-containing compounds in a single soil metagenome sample. c) Trend analysis results showing compounds with increasing (aqua) or decreasing (red) CPP with restoration. Clusters of points indicate energetically and chemically similar compounds. d) PCoA ordination showing CPP compositional differences (beta diversity)._

**Bioinformatic steps**

Metagenomics data were processed in several steps: (i) raw sequences were accessed/downloaded, (ii) QA/QC: sequences were inspected using FastQC and trimmed using Fastp, (iii) human genome sequences (GRCh38.p14/hg38) were removed from human gut datasets (in supplementary data, mouse genome sequences GRCm39 were removed from mouse gut datasets), (iv) functional profiles were derived using [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS), then (v) functional relative abundances were converted to CPP values. Metagenomic data processing to SUPER-FOCUS outputs, and subsequent conversion to CPP values were performed on [Pawsey Setonix](https://pawsey.org.au/systems/setonix/) linux high performance computers. Preparatory and intermediate steps, visualisation and data analysis were run on a local machine via [R code - main analyses](R-code-Part1-cpp_t2d-Main-analyses.R). We used the R [phyloseq](https://joey711.github.io/phyloseq/index.html) package for managing microbiome data. Supplementary R code analyses were performed for rarefying scenarios in the T2D datasets ([R code - rarefying analyses](R-code-Part2-cpp_t2d-Check-sequences-Rarefying-analyses-T2D-datasets.R)), and for preparation of additional preliminary validation datasets ([R code - supplementary analyses](R-code-Part3-cpp_t2d-Supplementary-data-analyses.R)). Bioinformatic steps run for each dataset are contained in folders in this repository. Note that folder/filepath structures used will need to be adjusted to run on other HPCs.  

**Software used**

Software versions used were: Python (>= v3.8.5), sratoolkit (v3.2.1), MG-RAST-Tools (v3.6.1), FastQC (v0.12.1), Fastp (v0.23.2), Bowtie2 (v2.4.1), Diamond (v0.9.19), SUPER-FOCUS (v0.34), seqtk (v1.5), R (HPC: 4.4.1, Local: v4.2.2).
