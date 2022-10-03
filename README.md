Rapidly Mixing Multiple-try Metropolis Algorithms for Model Selection
Problems
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# rapidMTM

<!-- badges: start -->
<!-- badges: end -->

This repository contains implementations of the paper [Rapidly Mixing
Multiple-try Metropolis Algorithms for Model Selection
Problems](https://arxiv.org/abs/2207.00689) by Hyunwoong Chang, Changwoo
Lee, Zhao Tang Luo, Huiyan Sang, and Quan Zhou, which is accepted at
[NeurIPS 2022](https://nips.cc/Conferences/2022).

We study multiple-try Metropolis (MTM) algorithm, which is an extension
of the Metropolis-Hastings (MH) algorithm by selecting the proposed
state among multiple trials according to some weight function
$w(y\,|\,x)$.

![mtm illustration](fig/mtm.PNG)

We show that multiple-try Metropolis (MTM) algorithm can achieve a
mixing time bound smaller than that of Metropolis-Hastings (MH)
algorithm by a factor of the number of trials under a general setting
applicable to high-dimensional model selection problems.

![Bayesian variable selection example](fig/fig3.PNG)

## Installation

You can install the development version of rapidMTM like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
