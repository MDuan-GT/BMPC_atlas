# BMPC_atlas
## Overview

This repository contains the analysis scripts to reproduce the main figures presented in the publication **Duan, Meixue et al. Understanding heterogeneity of human bone marrow plasma cell maturation and survival pathways by single-cell analyses. Cell Reports, Volume 42, Issue 7, 112682, [https://doi.org/10.1016/j.celrep.2023.112682](https://doi.org/10.1016/j.celrep.2023.112682)**

**Important Note on Data:** This analysis uses data from two sequencing runs. The second run was performed at a deeper sequencing depth aiming to improve the capture of non-immunoglobulin transcripts. However, this increased depth resulted in a sample index hopping issue, which was resolved in collaboration with 10x Genomics Support. The provided list of cleaned cell barcodes is in the data/clean_barcode folder.

## Processing Instructions
A critical pre-processing step is to filter the raw gene-cell count matrix using the provided cleaned barcodes from data/clean_barcode folder with the file name containing CLS055 and CLS058.

**For a practical implementation of this filtering step, please see the script for `Figure 1B`.**

