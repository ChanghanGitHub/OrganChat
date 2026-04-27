
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

``` r
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/github_Overview.png",
  out.width = "100%"
)
```

OrganChatDB integrates LS-receptor signaling, receptor-SE interactions,
and SE-target regulations for both human and mouse, as well as a
dictionary linking HMDB IDs to metabolite synonyms. In total, the
database catalogs 777 and 801 distinct LS molecules – encompassing
metabolites, peptide and nonpeptide hormones, and cytokines – for human
and mouse, respectively (Fig. 1a-b). The standard OrganChat workflow
consists of three key modules: (1) Data processing module; (2) OOC
inference module; (3) Multi-scale analysis module.

## Tutorial
