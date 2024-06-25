# Model per-gene variance in expression

![Unit tests](https://github.com/libscran/model_gene_variances/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/model_gene_variances/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/model_gene_variances/graph/badge.svg?token=JWV0I4WJX2)](https://codecov.io/gh/libscran/model_gene_variances)

## Overview

This repository contains functions to model the per-gene expression from a gene-by-cell matrix of (log-transformed) expression values.
Genes with high variance are considered to be more interesting and are prioritized for further analyses.
The code itself was originally derived from the [**scran** R package](https://bioconductor.org/packages/scran),
factored out into a separate C++ library for easier re-use.

## Quick start

Given a [`tatami::Matrix`](https://github.com/tatami-inc/tatami), the `model_gene_variances::compute()` function will compute several variance-related statistics for each gene:

```cpp
#include "scran/model_gene_variances.hpp"

tatami::Matrix<double, int>* ptr = some_data_source();

scran::model_gene_variances::Options opt;
auto res = scran::model_gene_variances::compute(ptr, opt);

res.means; // vector of means across genes.
res.variances; // vector of variances across genes.
res.fitted; // vector of fitted values of the mean-variance trend for each gene.
res.residuals; // vector of residuals from the trend.
```

Typically, the residuals are used for feature selection, as these account for non-trivial mean-variance trends in transformed count data.

```cpp
scran::choose_highly_variable_genes::Options copt;
copt.top = 5000;
auto chosen = scran::choose_highly_variable_genes::compute(
    res.residuals.size(), 
    res.residuals.data(), 
    copt
);
```

Users can also fit a trend directly to their own statistics.

```cpp
scran::fit_variance_trend::Options fopt;
fopt.span = 0.5;
fopt.minimum_mean = 1;
auto fit = scran::fit_variance_trend::compute(100, means, variances, fopt);
fit.fitted; // fitted values for all genes.
fit.residuals; // residuals values for all genes.
```

Check out the [reference documentation](https://libscran.github.io/model_gene_variances) for more details.

## Building projects

This repository is part of the broader [**libscran**](https://github.com/libscran/libscran) library,
so users are recommended to use the latter in their projects.
**libscran** developers should just use CMake with `FetchContent`:

```cmake
include(FetchContent)

FetchContent_Declare(
  scran_model_gene_variances 
  GIT_REPOSITORY https://github.com/libscran/model_gene_variances
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(scran_model_gene_variances)

# For executables:
target_link_libraries(myexe scran_model_gene_variances)

# For libaries
target_link_libraries(mylib INTERFACE scran_model_gene_variances)
```
