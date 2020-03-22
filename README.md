# EfficientFFBS
Efficient Forward Filtering and Backward Smoothing

Description: This repository contains implementations of the conventional and efficient forawrd filtering and backward smoothing.

Copyright (c) 2020 Behrad Soleimani All Rights Reserved

Contact: behrad@umd.edu

Citation: If you find these piece of codes helpful in your reserach, please cite the following technical report

-****

Date: March 5, 2020

Requirements: implemented in Matlab R2019a version, but should run on most versions.

Contens: 
> main.m:       Master script. 

> EFBS.m:       Efficient forward filtering and backward smoothing function.

> Filtering.m:  Conventional forward filtering and backward smoothing function.

> EfficientFFBS: Derivation and details of the algorithm.

Instructions: Simple and easy. Download all the codes in a directory and run main.m, that will generate one example described below. To use the functions individually, please look at the function descriptions.

Example:




| ![](Figs/Comparison.png) | 
|:--:| 
| Comparison of the conventional filtering and EfficientFFBS |



| ![](Figs/Quantile.png) | 
|:--:| 
| The estimation of EfficientFFBS with 95% quantile |
