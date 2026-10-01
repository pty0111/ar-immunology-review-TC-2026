# ar-immunology-review-TC-2026

This repository contains code for analysis in TC review article in Annual Review of Immunology.

## Overview

The repo contains code for analyzing and integrating scRNA-seq datasets of TCs.

## Setup

Code in this repo requires R 4.2.3 and Python 3.8 as well as the following dependencies:
- Seurat v4.4.0
- dplyr
- ggplot2

## Data availability

The datasets used in this repo are available at the following locations:

| Dataset Name             | Data accession                        | Publication                                                           |
|------------------------|------------------------|------------------------|
| TC_all_LN                | GSE294005                             | <https://www.science.org/doi/10.1126/science.adp0535>                 |
| Colonna                  | GSE312479                             | <https://www.science.org/doi/10.1126/sciimmunol.aed7439>              |
| TCp                    | GSE316304                             | <https://www.nature.com/articles/s41586-026-10198-z>                  |
| Sawa                     | GSE278161                             | <https://www.cell.com/immunity/fulltext/S1074-7613(26)00054-3>              |

The results of integration including Seurat R object are available at [Box](https://mskcc.box.com/s/sv0hs53icqy01l5ftwy8wdkurihg2eb1).
