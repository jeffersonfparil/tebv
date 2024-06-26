# tebv

Trial-estimated breeding values: best linear unbiased predictors of entries in single and multi-environment breeding trials

|**Build Status**|**License**|
|:--------------:|:---------:|
| <a href="https://github.com/jeffersonfparil/tebv/actions"><img src="https://github.com/jeffersonfparil/tebv/actions/workflows/r.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

A streamlined interface to calculate the breeding values of entries in breeding trials.

## Installation

```R
devtools::install_github("jeffersonfparil/tebv")
```

## Approach and principle

At the moment the approach will be to create a more or less homogenous interface for univariate analysis with plans to extend it to multivariate analysis.

We fit multiple *plausible* linear models via:

- canonical Henderson's equations (efficient if $n \ge p$) and,
- Newton-Raphson transformations (more efficient if $n \lt p$).
We then select the best fitting model based on which one on average has higher log-likelihood, lower AIC and lower BIC.

We fit 3 basic experimental design models:

- completely randomised design (CRD),
- randomised complete/incomplete block design (RBD), and
- spatial design with row and column effects (SPAT).

These models are implemented for:

- single environment trials, and
- multiple environment trials.

We use open-source tools as much as possible. This means we prefer `sommer` over `ASReml`. I would like to believe that being open-source where many other people can review and correct your code is a much better way of having confidence on the tool than having a lot of previous studies published using the tool.

## Naming convention

I like types, hence R being not really strongly typed I like to prefix my variable names with the type I wish that variable to remain to be throughout its lifetime:

- Generally no prefix for scalars, except for indexes (i.e., `idx_`) and booleans (i.e., `bool_`)
- `vec_` for vectors (a single lowercase letter can be used especially in the context of linear algebra, e.g. $y$, $\beta$ amd $\epsilon$ in $y = X\beta + \epsilon$)
- `mat_` for matrices (a single uppercase letter can be used especially in the context of linear algebra, e.g. $X$ in $\left( X'X \right)^{-1}$)
- `arr_` for arrays with more than 2 dimensions
- `list_` for lists
- `df_` for data frames
- `fn_` for functions
- `mod_` for models

The following are the names of the input data frame and its fields:

- The main input data frame per function is always called `df`.
- In univariate models, the trait is always called `y` with the corresponding function parameter `trait`, e.g. `trait="y"`.
- Environmental factor which may refer to different trials across time and space is always called `environ` with the corresponding function parameter `env`, e.g. `env="environ"`.
- Complete or incomplete blocks factor is always called `block` with the corresponding function parameter `block`, e.g. `block="block"`.
- Row factor is always called `row` with the corresponding function parameter `row`, e.g. `row="row"`.
- Column factor is always called `col` with the corresponding function parameter `col`, e.g. `col="col"`.

Model definitions:

- There should always be an explicit `1` to indicate the intecept in the model.
- The residual variance-covariance matrix should always be explicit even when inpedently and identically distributed, e.g. `rcov(units)`.


## Documentation

Write documentation or comments into your code a much as you can especially on top of your function definitions.

```R
devtools::document()
```

### Unit tests

Value modularity and write tests for each function definition.

See the main tests function and each Rscript for individual unit tests:

```R
devtools::check()
```

Or per module as:

```R
devtools::load_all()
library(sommer)
library(testthat)
source("tests/testthat/test-helpers.R")
source("tests/testthat/test-univariate_gx1.R")
source("tests/testthat/test-univariate_gxe.R")
```

Test new models in [`tests/testthat/test-helpers.R`](tests/testthat/test-helpers.R) and adjust parsing accordingly via: [`R/helpers.R::fn_henderson_vs_newtonraphson_fit()`](R/helpers.R).
