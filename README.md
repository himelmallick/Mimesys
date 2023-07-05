# Mimesys (Simulating Multi-omics Datasets with Arbitrary Association Structure)

## Introduction

MimESys provides a set of validated simulation modules to generate realistic multi-omics datasets facilitating targeted benchmarking of multi-omics analysis methods. 

## Installation

`MimESys` can be loaded using the following command (execute from within a fresh R session):
```r
install.packages('devtools')
library(devtools)
devtools::install_github("himelmallick/MimESys")
library(MimESys)
```

## Basic Usage

```r
MimESys(n_sample, n_feature, template_data)
```


## Citation

If you use `MimESys ` in your work, please cite the following paper:

Mallick H et al. (2022+). MimESys: Simulating Multi-omics Datasets with Arbitrary Association Structure. In Submission.


