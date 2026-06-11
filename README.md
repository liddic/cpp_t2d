# cpp_t2d
Examining shared trends in compound-associated functional capacities of degraded ecosystem soil microbiomes and type 2 diabetes (T2D) gut microbiomes


**! SITE UNDER CONSTRUCTION !**

**Contains code for the article:**

*Shared potential metabolism trends in degraded soils and type 2 diabetes gut microbiomes*

by Craig Liddicoat, Bart A. Eijkelkamp, Timothy R. Cavagnaro, Jake M. Robinson, Kiri Joy Wallace, Andrew D. Barnes, Garth Harmsworth, Damien J. Keating, Robert A. Edwards and Martin F. Breed

Preprint server: [BIORXIV/2025/642605](https://doi.org/10.1101/2025.03.11.642605)

**About this repository:** This repository documents code used in the above study, and is not intended for ongoing code development. Sub-folders contain separate workflow documentation and scripts relevant to the case study datasets described below.

**Data sources:** comprise soil metagenome data from five ecosystem quality gradients and two cohorts of untreated T2D versus normal healthy subjects

*Soil datasets*
- post-mining ecosystem restoration (data published by [Sun and Badgley 2019](https://www.sciencedirect.com/science/article/abs/pii/S0038071719301385?via%3Dihub))
- disturbed versus natural soils (data from the [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/australian-microbiome) previously examined by [Liddicoat et al 2024](https://doi.org/10.1016/j.scitotenv.2024.173543)
- prairie restoration (data published by [Mason et al 2023](https://academic.oup.com/femsec/article/99/11/fiad120/7288652))
- plantation succession (data published by [Zuo et al 2026](https://link.springer.com/article/10.1007/s11104-025-08246-0))
- vegetation succession (data published by [Wang et al 2024](https://academic.oup.com/ismecommun/article/4/1/ycae128/7833432))

*T2D datasets*
- Chinese cohort (data published by [Forslund et al 2015](https://www.nature.com/articles/nature15766) and [Qin et al 2012](https://www.nature.com/articles/nature11450))
- Swedish cohort (data published by [Forslund et al 2015](https://www.nature.com/articles/nature15766) and [Karlsson et al 2013](https://www.nature.com/articles/nature12198))

Study datasets are available from NCBI [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/) accessions PRJEB1786, PRJNA422434, PRJNA1215775, PRJNA1215778, PRJNA1215780, PRJNA1215781, PRJNA1080685; [MG-RAST](https://www.mg-rast.org/) project mgp16379; [JGI IMG](https://img.jgi.doe.gov/cgi-bin/m/main.cgi) Study ID Gs0144357, and the [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/australian-microbiome).

*Supplementary datasets*

Additional datasets were analysed to provide further validation of the new compound processing potential (CPP) method:
- Mice gut metagenomes raised in soil from desert, grassland and forest (data published by [Liu et al 2021](https://doi.org/10.1080/19490976.2020.1830699); NCBI Sequence read archive, accession: PRJNA542998)
- Soil metagenomes at 0h, 8h, 24h, and 48h following glucose amendment (data published by [Chuckran et al 2020](https://journals.asm.org/doi/10.1128/mra.00895-20); NCBI Sequence Read Archive accessions PRJNA539715, PRJNA539712, PRJNA539720, PRJNA539713, PRJNA539717, PRJNA539718, PRJNA539719, PRJNA539721, PRJNA539722, PRJNA539723, PRJNA539714, PRJNA539711, PRJNA539716.)
- Single-species gut bacteria culture metagenome standards (data published by [Amos et al 2020](https://link.springer.com/article/10.1186/s40168-020-00856-3) for the National Institute for Biological Standards and Control; NCBI Sequence Read Archive accession PRJNA622674)

**Method overview:** This study uses microbial metagenome functional profiling at the resolution of individual compounds, as outlined in [cpp3d](https://github.com/liddic/cpp3d) and the [manuscript](https://doi.org/10.1101/2025.03.11.642605). Metagenomics data were processed in several steps: (i) raw sequences were accessed/downloaded, (ii) QA/QC: sequences were inspected using FastQC and trimmed using Fastp, (iii) human genome sequences (GRCh38.p14/hg38) were removed from human gut datasets (in supplementary data, mouse genome sequences GRCm39 were removed from mouse gut datasets), then (iv) functional profiles were derived using [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS). Metagenomics data were processed on [Pawsey Setonix](https://pawsey.org.au/systems/setonix/) linux high performance computers. SUPER-FOCUS results were downloaded to a local machine for further visualisation and analysis via [R code - TO BE UPDATED!](link-under-construction.R). These analyses represent a new approach to translate community-scale microbiome functional potential relative abundances down to potential metabolism at the scale of individual compounds. Bioinformatic steps run for each dataset are contained in folders in this repository. Note that folder/filepath structures used will need to be adjusted to run on other HPCs. Software versions used were: Python (v3.8.5), FastQC (v0.11.9), Fastp (v0.23.2), Bowtie2 (v2.4.1), Diamond (v0.9.19), SUPER-FOCUS (v0.0.0), R (v4.2.2).
