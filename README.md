# vcmasf
This package implements the one-step and two-step spline fitting with adaptive knots selection for varying coefficient modeling, as well as variable selection method for sparse high dimension problem. More details are in our paper `Varying Coefficient Modeling via Adaptive Spline Fitting`.



## Install the package 'vcmasf'

The package requires `Rcpp`, `RcppArmadillo`, `splines`, `pracma` before installation. To install the package from Github, please use

```
devtools::install_github('wangxf0106/vcmasf').
```

When installing package on Mac, for error messages like `ld: library not found for -lgfortran`, you can follow [RCPP, RCPPARMADILLO AND OS X MAVERICKS "-LGFORTRAN" AND "-LQUADMATH" ERROR](https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/) to install gfortran.

## Examples

The R code in `example/example1.R` shows the one-step and two-step fitting methods with adaptive knots selection.

The R code in `example/example2.R` presents the variable selection method for sparse high dimension problem.

The R code in `example/example3.R` applies the two-step fitting method to the COVID-19 case number and environmental factor dataset.

## Acknowledgement

We are grateful to the National Oceanic and Atmospheric Administration Regional Climate Centers, Northeast Regional Climate Center at Cornell University for kindly sharing the meteorological data. We are thankful to the Environmental Protection Agency and the Department of Health, New York State for publishing the air quality data and daily infected records online. 
