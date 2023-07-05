# MimeSys (Simulating Multi-omics Datasets with Arbitrary Association Structure)

## Introduction

MimeSys provides a set of validated simulation modules to generate realistic multi-omics datasets facilitating targeted benchmarking of multi-omics analysis methods. 

## Installation

`MimeSys` can be loaded using the following command (execute from within a fresh R session):
```r
install.packages('devtools')
library(devtools)
devtools::install_github("himelmallick/MimeSys")
library(MimeSys)
```

## Basic Usage

```r
MimeSys(n_sample, n_feature, template_data)
```


## Citation

If you use `MimeSys ` in your work, please cite the following paper:

Mallick H et al. (2023+). MimESys: Simulating Multi-omics Datasets with Arbitrary Association Structure. In Submission.


