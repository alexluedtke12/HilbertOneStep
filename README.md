<!--
---
title: "HilbertOneStep R Package"
output: rmarkdown::github_document
---
-->



# Description

The HilbertOneStep package implements one-step estimators of Hilbert-valued parameters, as evaluated in the simulation study from:

"One-Step Estimation of Differentiable Hilbert-Valued Parameters" by A. Luedtke and I. Cheung [[link](https://arxiv.org/abs/XXXXXXXX)]

# Installation

To install the package using [devtools](https://www.rstudio.com/products/rpackages/devtools/), run the following commands:


```r
library(devtools)
devtools::install_github("alexluedtke12/HilbertOneStep")
```

# Examples

This package includes four key functions. The first two estimate the counterfactual density under treatment A=1, assuming there are no unmeasured confounders and the positivity assumption holds. The first function implements a regularized one-step estimator that is consistent within a nonparametric model under conditions. See Luedtke and Chung (2023) for rate of convergence guarantees for this estimator.

```r
# load the package
library(HilbertOneStep)

# sample size
n = 500

# simulate data
dat = sim_data(n,setting='zeroBothSides')
W = dat$W
A = dat$A
Y = dat$Y

cv_out = cv_density(W,A,Y,ngrid=500,num_fold_final=2)
cv_out$which_best
ygrid = seq(0.0025,0.9975,by=0.0025)
plot(ygrid,cv_out$best_fun(ygrid),xlab='y',ylab='Counterfactual Density Estimate',type='l')
```
<br>

If it is known that the counterfactual density is bandlimited, a non-regularized one-step estimator of the counterfactual density is available. This estimator achieves a parametric rate of convergence, with mean integrated squared error converging to zero at a rate of 1/n.

```r
# sample size
n = 500

# simulate data
dat = sim_data_bandlimited(n)
W = dat$W
A = dat$A
Y = dat$Y

bl_out = bandlimited_density(W,A,Y,ngrid=500,num_fold=2,num_boot=2000,alpha=seq(0.01,0.99,by=0.01),b=2)
ygrid = seq(-15,15,by=0.01)
plot(ygrid,bl_out$onestep_fun(ygrid),xlab='y',ylab='Counterfactual Density Estimate',type='l')
```
<br>

The package also implements tests to determine whether the counterfactual distributions of the outcome under interventions that set A=1 and A=0 are the same.

```r
# sample size
n = 250

# simulate data from alternative
dat = sim_data(n,setting='nonzeroBothSides',cos_parameter=1)
W = dat$W
A = dat$A
Y = dat$Y

# a test based on a statistic that contrasts the counterfactual density functions under A=1 and A=0.
density_test(W,A,Y,num_fold=2)
# a test based on the maximum mean discrepancy between the counterfactual distributions under A=1 and A=0.
mmd_test(W,A,Y,num_fold=2)


# simulate data from null
dat = sim_data(n,setting='nonzeroBothSides',cos_parameter=0)
W = dat$W
A = dat$A
Y = dat$Y

# a test based on a statistic that contrasts the counterfactual density functions under A=1 and A=0.
density_test(W,A,Y,num_fold=2)
# a test based on the maximum mean discrepancy between the counterfactual distributions under A=1 and A=0.
mmd_test(W,A,Y,num_fold=2)
```
