
# Expression Analysis in tumors for Paper "Zbtb38 transcriptionally activates XIAP to regulate apoptosis in development and cancer"

This repository contains all scripts and materials for reproducing the analyses in the paper "Zbtb38 transcriptionally activates XIAP to regulate apoptosis in development and cancer".

## Directory Structure

    .
    ├── scripts/
    │   ├── data_download.r
    │   ├── merging_data.r
    │   ├── figureA.r
    │   ├── figureB.r
    │   ├── figureCDE.r
    │   ├── figureF.r
    │   ├── figureG.r
    │   ├── figureH.r
    │   └── Sfigures.r
    ├── data/
    │   ├── MERAV_data/
    │   └── TCGA_data/
    │       └── all_cancer_APOPTOTIC_correlation_withZBTB38_modified
    ├── figures/
    └── README.md

## Requirements

- R: 4.3.2
    - TCGAbiolinks: 2.30.4
    - here: 1.0.1
    - SummarizedExperiment: 1.32.0
    - dplyr: 1.1.4
    - DESeq2: 1.42.1
    - ggplot2: 3.5.0
    - scales: 1.3.0
    - irlba: 2.3.5.1
    - Rtsne: 0.17
    - ggrepel: 0.9.5
    - ggpubr: 0.6.0
    - readr: 2.1.5
    - viridis: 0.6.5


## How to Run

1. Download "TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes" from
UCSC Xena Browser(https://xenabrowser.net/datapages/) and place it in "data/TCGA_data".
2. Download "all_genes_batch_adjusted_RAW.txt" from Metabolic gEne RApid Visualizer(http://merav.wi.mit.edu/) and place it in "data/MERAV_data".
3. From the root directory, run each script in the following order using Rscript:
    ```sh
    Rscript scripts/data_download.r
    Rscript scripts/merging_data.r
    Rscript scripts/figureA.r
    Rscript scripts/figureB.r
    Rscript scripts/figureCDE.r
    Rscript scripts/figureF.r
    Rscript scripts/figureG.r
    Rscript scripts/figureH.r
    Rscript scripts/Sfigures.r
    ```
4. Figures will be saved in `figures/`.

## Script–Figure Correspondence

| Script/File        | Output File         | Corresponds to Figure in Paper  |
|--------------------|---------------------|---------------------------------|
| scripts/figureA.r  | figures/figureA.pdf | Figure 5A                       |
| scripts/figureB.r  | figures/figureB.pdf | Figure 5B                       |
| scripts/figureCDE.r| figures/figureC.pdf | Figure 5C                       |
| scripts/figureCDE.r| figures/figureD.pdf | Figure 5D                       |
| scripts/figureCDE.r| figures/figureE.pdf | Figure 5E                       |
| scripts/figureF.r  | figures/figureF.pdf | Figure 5F                       |
| scripts/figureG.r  | figures/figureG.pdf | Figure 5G                       |
| scripts/figureH.r  | figures/figureH.pdf | Figure 5H                       |
| scripts/Sfigures.r | figures/SFigure1.pdf| Supplementary Figures           |

## Note on Data Reproducibility

All TCGA data used in this study were downloaded via TCGAbiolinks (v2.30.4) on 2025-02-14.
Because the TCGA and GDC data portals are continuously updated, future downloads may not exactly match the dataset used here.

## Citation
If you use these scripts or data, please cite our paper as follows: 

Toshiaki Shigeoka, Hiroyuki Nagaoka, Nunuk Aries Nurulita, Shogo Tada, Yasumasa Bessho, Yasumasa Ishida, Eishou Matsuda, Zbtb38 transcriptionally activates XIAP to regulate apoptosis in development and cancer, in submission. 

## Contact
Toshiaki Shigeoka
Eishou Matsuda
Address for correspondence in the paper 

