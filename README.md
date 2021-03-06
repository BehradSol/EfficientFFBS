# EfficientFFBS
Efficient Forward Filtering and Backward Smoothing

Description: This repository contains implementations of the conventional and efficient forawrd filtering and backward smoothing.

Copyright (c) 2020 Behrad Soleimani All Rights Reserved

Contact: behrad@umd.edu

Citation: If you find these piece of codes helpful in your reserach, please cite the following paper

- B. Soleimani, P. Das, J. Kulasingham, J. Z. Simon and B. Babadi, "Granger Causal Inference from Indirect Low-Dimensional Measurements with Application to MEG Functional Connectivity Analysis," 2020 54th Annual Conference on Information Sciences and Systems (CISS), Princeton, NJ, USA, 2020, pp. 1-5, [doi: 10.1109/CISS48834.2020.1570617418](https://ieeexplore.ieee.org/abstract/document/9086218).

Date: March 5, 2020

Requirements: implemented in Matlab R2019a version, but should run on most versions.

Contents: 
> main.m:       **Master script**. 

> EFBS.m:       **Efficient forward filtering and backward smoothing function**.

> Filtering.m:  **Conventional forward filtering and backward smoothing function**.

> EfficientFFBS.pdf: **Derivation and details of the algorithm**.

Instructions: Simple and easy. Download all the codes in a directory and run main.m, that will generate one example described below. To use the functions individually, please look at the function descriptions. The derivations and details are also explained in .pdf file.

Example:

We consider the following observation model

<p align="center">
  <img src="https://user-images.githubusercontent.com/59627073/81014091-39efaa00-8e2a-11ea-8640-24d4fb30b3cd.jpg">
</p>
where x<sub>t</sub> and y<sub>t</sub> represent the observation and source vectors at t-th time sample, respectively. The underlying source dynamic is modeled via a vector autoregressive process, VAR(p), as
<p align="center">
  <img src="https://user-images.githubusercontent.com/59627073/81014376-bd110000-8e2a-11ea-91e4-e41cb0ac6543.jpg">
</p>
where e<sub>t</sub> shows the (external) stimuli vector corresponding to the t-th time sample. The goal is to obtain the non-causal beliefs, i.e. 

<p align="center">
  <img src="https://user-images.githubusercontent.com/59627073/81014834-969f9480-8e2b-11ea-975c-b01266ac2f14.jpg">
</p>

In this example, we assume that there are N<sub>y</sub>=3 observations, N<sub>x</sub>=5 sources, and T=50 time samples. The underlying source dynamic is considered as VAR(1). 

In Fig.1, the estimated version of the source #1, i.e. non-causal beliefs p(x<sub>t</sub> | y<sub>1:T</sub>), and ground truth are compared. As can be seen, the performance of the conventional filtering scheme (Kalman filtering) is almost the same as the EfficientFFBS. 


| ![](Figs/Comparison.png) | 
|:--:| 
| Fig 1. Comparison of the conventional filtering and EfficientFFBS |

In Fig.2, the 95% quantile of the estimation (via EfficientFFBS) is depicted to demonstrate the accuracy of the estimation.

| ![](Figs/Quantile.png) | 
|:--:| 
| Fig 2. The estimation of EfficientFFBS with 95% quantile |
