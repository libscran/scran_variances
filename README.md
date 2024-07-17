# Model per-gene variance in expression

![Unit tests](https://github.com/libscran/scran_variances/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/scran_variances/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/scran_variances/graph/badge.svg?token=JWV0I4WJX2)](https://codecov.io/gh/libscran/scran_variances)

## Overview

This repository contains functions to model the per-gene expression from a gene-by-cell matrix of (log-transformed) expression values.
Genes with high variance are considered to be more interesting and are prioritized for further analyses.
The code itself was originally derived from the [**scran** R package](https://bioconductor.org/packages/scran),
factored out into a separate C++ library for easier re-use.

## Quick start

Given a [`tatami::Matrix`](https://github.com/tatami-inc/tatami) of log-expression values for each gene in each cell,
we can compute the per-gene variances and model the trend with respect to the mean across genes:

```cpp
#include "scran_variances/scran_variances.hpp"

std::shared_ptr<tatami::Matrix<double, int> > mat = some_data_source();

scran_variances::ModelGeneVariancesOptions opt;
auto res = scran_variances::model_gene_variances(*mat, opt);

res.means; // vector of means across genes.
res.variances; // vector of variances across genes.
res.fitted; // vector of fitted values of the mean-variance trend for each gene.
res.residuals; // vector of residuals from the trend.
```

Typically, the residuals are used for feature selection, as these account for non-trivial mean-variance trends in transformed count data.

```cpp
scran_variances::ChooseHighlyVariableGenesOptions copt;
copt.top = 5000;
auto chosen = scran_variances::choose_highly_variable_genes_index(
    res.residuals.size(), 
    res.residuals.data(), 
    copt
);

// Create the HVG submatrix for downstream analysis.
auto hvg_subset = tatami::make_DelayedSubset(mat, chosen, /* by_row = */ true);
```

Users can also fit a trend directly to their own statistics.

```cpp
scran_variances::FitVarianceTrendOptions fopt;
fopt.span = 0.5;
fopt.minimum_mean = 1;
auto fit = scran_variances::fit_variance_trend(100, means, variances, fopt);
fit.fitted; // fitted values for all genes.
fit.residuals; // residuals values for all genes.
```

Check out the [reference documentation](https://libscran.github.io/scran_variances) for more details.

## Building projects

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  scran_variances
  GIT_REPOSITORY https://github.com/libscran/scran_variances
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(scran_variances)
```

Then you can link to **scran_variances** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe libscran::scran_variances)

# For libaries
target_link_libraries(mylib INTERFACE libscran::scran_variances)
```

### CMake with `find_package()`

```cmake
find_package(libscran_scran_variances CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE libscran::scran_variances)
```

To install the library, use:

```sh
mkdir build && cd build
cmake .. -DSCRAN_VARIANCES_TESTS=OFF
cmake --build . --target install
```

By default, this will use `FetchContent` to fetch all external dependencies.
If you want to install them manually, use `-DSCRAN_VARIANCES_FETCH_EXTERN=OFF`.
See the tags in [`extern/CMakeLists.txt`](extern/CMakeLists.txt) to find compatible versions of each dependency.

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt), which also need to be made available during compilation.
