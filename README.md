# tebv
Trial-estimated breeding values: best linear unbiased predictors of entries in single and multi-environment breeding trials

|**Build Status**|**License**|
|:--------------:|:---------:|
| <a href="https://github.com/jeffersonfparil/tebv/actions"><img src="https://github.com/jeffersonfparil/tebv/actions/workflows/r.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

A streamlined interface to calculate the breeding values of entries in breeding trials.

## Approach and principle

At the moment the approach will be to create a more or less homogenous interface for univariate analysis with plans to extend it to multivariate analysis.

We fit multiple *plausible* linear models via:
- canonical Henderson's equations (efficient if $n \ge p$) and,
- Newton-Raphson transformations (more efficient if $n \lt p$).
We then select the best fitting model based on which one on average has higher log-likelihood, lower AIC and lower BIC.

We fit 3 basic experimental design models:
- completely randomised design (CRD),
- randomised complete block design (RCBD), and
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
- `df_` for data.frames
- `fn_` for functions
- `mod_` for models

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
