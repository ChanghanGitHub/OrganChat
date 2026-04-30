
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OrganChat: inferring organ-organ communication through single-cell data

OrganChat is a computational method that uses single-cell multi-modal
data together with a newly constructed long-range signal (LS)-mediated
communication database (OrganChatDB) to infer organ-organ
communications.

## Installation

You can install the development version of OrganChat from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ChanghanGitHub/OrganChat")
```

## Overview

<img width="536" alt="image" src="man/figures/github_Overview.png"/>

OrganChatDB integrates LS-receptor signaling, receptor-SE interactions,
and SE-target regulations for both human and mouse, as well as a
dictionary linking HMDB IDs to metabolite synonyms. In total, the
database catalogs 777 and 801 distinct LS molecules – encompassing
metabolites, peptide and nonpeptide hormones, and cytokines – for human
and mouse, respectively (Fig. 1a-b). The standard OrganChat workflow
consists of three key modules: (1) Data processing module; (2) OOC
inference module; (3) Multi-scale analysis module.

All processed data and some saved OrganChat objects can be download
from: <https://doi.org/10.5281/zenodo.19042679>.

## Tutorials

[Inferring, analyzing, and comparing the organ-organ communication (OOC)
between different conditions using
OrganChat.](https://htmlpreview.github.io/?https://github.com/ChanghanGitHub/OrganChat/blob/master/vignettes/Comparison_Analysis_Tutorial.html)

## Examples

### Benchmarking and perturbation-based validation using a spatial multi-omics human gastric cancer dataset: [link to the folder](Example/TumorMeta_Sun_2023)

1)  [Organ-organ communicarion analysis between two
    “pseudo-organs”.](Example/TumorMeta_Sun_2023/TumorMeta_Sun_2023_OC.Rmd)

2)  [Benchmarking analysis against CellChat, ICELLNET, and
    exFINDER](Example/TumorMeta_Sun_2023/TumorMeta_Sun_2023_Benchmarking.Rmd)

3)  [Perturbation-based
    validation.](Example/TumorMeta_Sun_2023/TumorMeta_Sun_2023_Perturbation_CAOC.Rmd)

### Application to a multi-organ insulin-resistance and aging-related system: [link to the folder](Example/Insulin_UCI_2025_mouse)

1)  [Data processing and metabolite flux
    inference.](Example/Insulin_UCI_2025_mouse/MultiOrgan_Mouse_UCI_DataProcessingMetaboliteAnalysis.Rmd)

2)  [Evaluating OrganChat predictions using a two-organ
    insulin-resistance
    system.](Example/Insulin_UCI_2025_mouse/MultiOrgan_Mouse_UCI_Validation_Liver_Adipose_OS_OR.Rmd)

3)  [Further evaluation of an extended three-organ system for
    aging-related
    diabetes.](Example/Insulin_UCI_2025_mouse/MultiOrgan_Mouse_UCI_MOOC_T2D.Rmd)

### Application to a five-organ communication system for immune-mediated diseases: [link to the folder](Example/Sandra_CellRepMed_2023)

1)  [Processing the original data and inferring metabolite
    flux.](Example/Sandra_CellRepMed_2023/Sandra_CellRepMed_2023_DataProcessing_MetaboliteAnalysis.Rmd)

2)  [Comparison analysis via OrganChat between conditions for: joint and
    muscle; joint and spleen; muscle and spleen; skin and
    muscle.](Example/Sandra_CellRepMed_2023/Sandra_CellRepMed_2023_CAOC.Rmd)

3)  [Multi-organ OrganChat (MOOC) Analysis of OrganChat between
    conditions for: joint, muscle, lung, skin,
    spleen.](Example/Sandra_CellRepMed_2023/Sandra_CellRepMed_2023_MOOC.Rmd)

### OrganChat reveals cross-organ communication for human multiple-disease systems: [link to the folder](Example/HumanAtlas_Science_2022)

1)  [Processing the original
    data.](Example/HumanAtlas_Science_2022/HumanAtlas_Science_2022_DataProcessing.Rmd)

2)  [Inferring metabolite flux
    data.](Example/HumanAtlas_Science_2022/HumanAtlas_Science_2022_MetaboliteAnalysis.Rmd)

3)  [Multi-organ OrganChat analysis between spleen, small intestine,
    lymph node, lung, and
    vasculature.](Example/HumanAtlas_Science_2022/HumanAtlas_Science_2022_MultiOrganAnalysis.Rmd)

4)  [Comparison analysis between two donors for: small intestine and
    vasculature, lung and vasculature, lymph node and vasculature,
    respectively.](Example/HumanAtlas_Science_2022/HumanAtlas_Science_2022_ComparisonAnalysis.Rmd)
