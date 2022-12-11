# vcmasf
This package implements the equidistant spline fitting, global and predictor-specific spline fitting with adaptive knots selection for varying coefficient modeling, as well as variable selection method for sparse high dimension problem. More details are in our paper [1].

## Install the package 'vcmasf'

The package requires `Rcpp`, `RcppArmadillo`, `splines`, `pracma` before installation. To install the package from Github, please use

```
devtools::install_github('wangxf0106/vcmasf').
```

When installing package on Mac, for error messages like `ld: library not found for -lgfortran`, you can follow [this guidance](https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/) to install gfortran.

## Examples
We provide example codes to replicate the simulation studies and real data applications in [1].

The R code in `example/example1.R` presents the simulation study in Section 4.1: Simulation study for adaptive spline fitting.

The R code in `example/example2.R` presents the simulation study in Section 4.2: Simulation study for variable selection.

The R code in `example/example3.R` presents the real data application in Section 5.1: Environmental factors and COVID-19.

The R code in `example/example4.R` presents the real data application in Section 5.2: Boston housing data.

## Acknowledgement

We are grateful to the National Oceanic and Atmospheric Administration Regional Climate Centers, Northeast Regional Climate Center at Cornell University for kindly sharing the meteorological data. We are thankful to the Environmental Protection Agency and the Department of Health, New York State for publishing the air quality data and daily infected records online. 

## References
<a id="1"> [1]</a> Wang, X., Jiang, B., Liu, J. S., Varying Coefficient Model via Adaptive Spline Fitting, [2022](https://arxiv.org/abs/2201.10063).
