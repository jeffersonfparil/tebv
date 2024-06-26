% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_gxe.R
\name{fn_assess_GXE}
\alias{fn_assess_GXE}
\title{Assess genotype-by-environment (GXE) interactions by calculating pairwise trait correlations between environments.
The trait values are simply aggregated using arithmetic mean and the pairwise correlations are calculated.}
\usage{
fn_assess_GXE(df, trait = "y", id = "gen", env = "env", verbose = FALSE)
}
\arguments{
\item{df}{data frame containing the model variables}

\item{trait}{name of the continuous numeric response variable in the data frame}

\item{id}{name of the entry field in the data frame}

\item{env}{name of the environment field in the data frame}

\item{verbose}{show scatterplots?}
}
\value{
mat_corr: correlations (Pearson's) between environments
mat_pval: p-values of the correlations
}
\description{
Assess genotype-by-environment (GXE) interactions by calculating pairwise trait correlations between environments.
The trait values are simply aggregated using arithmetic mean and the pairwise correlations are calculated.
}
\examples{
df = fn_simulate_gx1(design="crd")
out = fn_assess_GXE(df=df, trait="y", id="gen", env="env", verbose=TRUE)
}
