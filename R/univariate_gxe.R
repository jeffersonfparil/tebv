#' Assess genotype-by-environment (GXE) interactions by calculating pairwise trait correlations between environments.
#' The trait values are simply aggregated using arithmetic mean and the pairwise correlations are calculated.
#' 
#' @param df data frame containing the model variables
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field in the data frame
#' @param env name of the environment field in the data frame
#' @param verbose show scatterplots?
#' @returns
#' mat_corr: correlations (Pearson's) between environments
#' mat_pval: p-values of the correlations
#' @examples
#' df = fn_simulate_gx1(design="crd")
#' out = fn_assess_GXE(df=df, trait="y", id="gen", env="env", verbose=TRUE)
#' @export
fn_assess_GXE = function(df, trait="y", id="gen", env="env", verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # df = fn_simulate_gxe(design="crd")
    # trait="y"; id="gen"; env="env"; verbose=TRUE
    ### TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", env, "`"))))) {
        print(paste0("Error: env: ", env, " does not exist in the input dataframe."))
        return(NA)
    }
    y = eval(parse(text=paste0("df$`", trait, "`")))
    gen = eval(parse(text=paste0("df$`", id, "`")))
    environ = eval(parse(text=paste0("df$`", env, "`")))
    if (sum(!is.na(y))==0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    if (sum(!is.na(gen))==0) {
        print(paste0("Error: all genotype IDs are missing for the ID field: ", id, "."))
        return(NA)
    }
    if (sum(!is.na(environ))==0) {
        print(paste0("Error: all environment IDs are missing for the environment field: ", env, "."))
        return(NA)
    }
    df = data.frame(y=y, gen=gen, environ=environ)
    df_mu = stats::aggregate(y ~ environ + gen, FUN=mean, na.rm=TRUE, data=df)
    df_sd = stats::aggregate(y ~ environ + gen, FUN=stats::sd, na.rm=TRUE, data=df)
    vec_gen = unique(df_mu$gen)
    vec_environ = unique(df_mu$environ)
    n = length(vec_gen)
    m = length(vec_environ)
    mat_corr = matrix(1, nrow=m, ncol=m); rownames(mat_corr) = colnames(mat_corr) = vec_environ
    mat_pval = matrix(0, nrow=m, ncol=m); rownames(mat_pval) = colnames(mat_pval) = vec_environ
    for (i in 1:(m-1)) {
        for (j in (i+1):m) {
            # i = 1; j = 2
            idx_i = which(df_mu$environ == vec_environ[i])
            idx_j = which(df_mu$environ == vec_environ[j])
            y_i = df_mu$y[idx_i]
            y_j = df_mu$y[idx_j]
            corr_test = stats::cor.test(y_i, y_j)
            mat_corr[i, j] = mat_corr[j, i] = corr_test$estimate
            mat_pval[i, j] = mat_pval[j, i] = corr_test$p.value
            if (verbose) {
                print("###################################")
                print(paste0(vec_environ[i], " X ", vec_environ[j]))
                print(paste0("corr=", round(corr_test$estimate*100), "% (pval=", round(corr_test$p.value, 4), ")"))
                txtplot::txtplot(y_i, y_j)
            }
        }
    }
    return(list(mat_corr=mat_corr, mat_pval=mat_pval))
}

#' Fit a GXE model to extract the best linear unbiased predictors for the genotype values,
#' assuming a completely randomised design (CRD) in multiple environments.
#' We use this design if per environment, we do not expect heterogeneity in each trial area, e.g. small trial and highly-controlled conditions.
#' 
#' @param df data frame containing the model variables
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field in the data frame
#' @param env name of the environment field in the data frame
#' @param tolParInv diagonal loading of the variance-covariance matrix prior to inversion
#' @param skip_algo_based_on_dimensions if there are potentially more effects to be estimated than the number of observations then skip Henderson's fitting in favour of the Newton-Raphson transformations, and vice-versa
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values per environment (fields: id, env, BLUPs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the random effects
#' @examples
#' df = fn_simulate_gx1(design="crd")
#' out = fn_GXE_CRD_BLUPs(df=df, trait="y", id="gen", env="env", tolParInv=0.01, 
#'      skip_algo_based_on_dimensions=TRUE, verbose=TRUE)
#' @export
fn_GXE_CRD_BLUPs = function(df, trait="y", id="gen", env="env", tolParInv=0.01, skip_algo_based_on_dimensions=TRUE, verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gxe(design="crd")
    # trait="y"; id="gen"; env="env"; tolParInv=0.01; skip_algo_based_on_dimensions=TRUE; verbose=TRUE
    ### TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", env, "`"))))) {
        print(paste0("Error: env: ", env, " does not exist in the input dataframe."))
        return(NA)
    }
    ### Create consistently named response and explanatory variables for code readability below
    df_gxe_crd = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = eval(parse(text=paste0("df$", id))), 
        environ = eval(parse(text=paste0("df$", env)))
    )
    if (sum(!is.na(df_gxe_crd$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    mat_design_matrix_crd_gxe = stats::model.matrix(y ~ id:environ, data=df_gxe_crd)
    if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_crd_gxe) < ncol(mat_design_matrix_crd_gxe))) {
        ### if n < p (more coefficients to predict) then use Newton-Raphson transformations via mmer()
        skip_newtonrap = FALSE
        skip_henderson = TRUE
        mod_henderson_1 = mod_henderson_2 = mod_henderson_3 = mod_henderson_4 = mod_henderson_5 = mod_henderson_6 = NA
    } else if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_crd_gxe) >= ncol(mat_design_matrix_crd_gxe))) {
        ### if n >= p (less coefficients to predict) then use Henderson's mixed model equations via mmec()
        skip_newtonrap = TRUE
        skip_henderson = FALSE
        mod_newtonrap_1 = mod_newtonrap_2 = mod_newtonrap_3 = mod_newtonrap_4 = mod_newtonrap_5 = mod_newtonrap_6 = NA
    } else {
        skip_newtonrap = FALSE
        skip_henderson = FALSE
    }
    if (!skip_newtonrap) {
        if (verbose) {
            print("Fitting models via Newton-Raphson algorithm")
        }
        ###m Newton-Raphson models
        mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ environ:id, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + environ:id, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ vsr(usr(environ), id), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=2)), id), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_6 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=1)), id), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    if (!skip_henderson) {
        if (verbose) {
            print("Fitting models via Henderson's mixed model equations")
        }
        ### Henderson's models
        mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ environ:id, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + environ:id, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ vsc(usc(environ), isc(id)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=2)), isc(id)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_6 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=1)), isc(id)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    ### Identify the best model
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2, mod_henderson_3, mod_henderson_4, mod_henderson_5, mod_henderson_6), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2, mod_newtonrap_3, mod_newtonrap_4, mod_newtonrap_5, mod_newtonrap_6),
        verbose=verbose)
    list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=TRUE, verbose=verbose)
    if (is.na(list_u_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    if (verbose) {
        print(paste0("'Best fit': {", 
            "'algorithm': '", list_u_V_fitstats$algorithm, "', ", 
            "'model': {", 
                paste(paste0("'", names(list_u_V_fitstats$model), "': '", list_u_V_fitstats$model, "', "), collapse=" "),
            "}",
        "}"))
    }
    ### Extract breeding values per environment
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_6, mod_newtonrap=NA, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod_newtonrap_6, verbose=verbose)
    df_BLUPs_GXE = fn_extract_gxe_breeding_values(list_u_V_fitstats)
    return(list(
        df_effects=df_BLUPs_GXE, 
        loglik=list_u_V_fitstats$loglik,
        AIC=list_u_V_fitstats$AIC,
        BIC=list_u_V_fitstats$BIC,
        algorithm=list_u_V_fitstats$algorithm,
        model=list_u_V_fitstats$model,
        V=list_u_V_fitstats$V))
}

#' Fit a GXE model to extract the best linear unbiased predictors for the genotype values,
#' assuming a randomised complete block design (RCBD; each block is a complete replication of all genotypes) in multiple environments.
#' We use this design when we expect heterogeneity along a single direction in each trial area, e.g. field trial along a slope.
#' 
#' @param df data frame containing the model variables
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field in the data frame
#' @param env name of the environment field in the data frame
#' @param block name of the complete or incomplete block field in the data frame
#' @param tolParInv diagonal loading of the variance-covariance matrix prior to inversion
#' @param skip_algo_based_on_dimensions if there are potentially more effects to be estimated than the number of observations then skip Henderson's fitting in favour of the Newton-Raphson transformations, and vice-versa
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values per environment (fields: id, env, BLUPs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the random effects
#' @examples
#' df = fn_simulate_gx1(design="rbd")
#' out = fn_GXE_RBD_BLUPs(df=df, trait="y", id="gen", env="env", block="rep", tolParInv=0.01, 
#'      skip_algo_based_on_dimensions=TRUE, verbose=TRUE)
#' @export
fn_GXE_RBD_BLUPs = function(df, trait="y", id="gen", env="env", block="rep", tolParInv=0.01, skip_algo_based_on_dimensions=TRUE, verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gxe(design="rbd")
    # trait="y"; id="gen"; env="env"; block="rep"; tolParInv=0.01; skip_algo_based_on_dimensions=TRUE; verbose=TRUE
    ### TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", env, "`"))))) {
        print(paste0("Error: env: ", env, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", block, "`"))))) {
        print(paste0("Error: block: ", block, " does not exist in the input dataframe."))
        return(NA)
    }
    ### Create consistently named response and explanatory variables for code readability below
    df_gxe_rbd = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = eval(parse(text=paste0("df$", id))), 
        environ = eval(parse(text=paste0("df$", env))),
        block = eval(parse(text=paste0("df$", block)))
    )
    if (sum(!is.na(df_gxe_rbd$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    mat_design_matrix_rbd_gxe = stats::model.matrix(y ~ id:environ:block, data=df_gxe_rbd)
    if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_rbd_gxe) < ncol(mat_design_matrix_rbd_gxe))) {
        ### if n < p (more coefficients to predict) then use Newton-Raphson transformations via mmer()
        skip_newtonrap = FALSE
        skip_henderson = TRUE
        mod_henderson_1 = mod_henderson_2 = mod_henderson_3 = mod_henderson_4 = mod_henderson_5 = mod_henderson_6 = NA
    } else if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_rbd_gxe) >= ncol(mat_design_matrix_rbd_gxe))) {
        ### if n >= p (less coefficients to predict) then use Henderson's mixed model equations via mmec()
        skip_newtonrap = TRUE
        skip_henderson = FALSE
        mod_newtonrap_1 = mod_newtonrap_2 = mod_newtonrap_3 = mod_newtonrap_4 = mod_newtonrap_5 = mod_newtonrap_6 = NA
    } else {
        skip_newtonrap = FALSE
        skip_henderson = FALSE
    }
    if (!skip_newtonrap) {
        if (verbose) {
            print("Fitting models via Newton-Raphson algorithm")
        }
        ###m Newton-Raphson models
        mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ environ:id + environ:block, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + environ:id + environ:block, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ vsr(usr(environ), id) + vsr(usr(environ), block), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id) + vsr(usr(environ), block), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=2)), id) + vsr(usr(environ), block), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_6 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=1)), id) + vsr(usr(environ), block), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    if (!skip_henderson) {
        if (verbose) {
            print("Fitting models via Henderson's mixed model equations")
        }
        ### Henderson's models
        mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ environ:id + environ:block, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + environ:id + environ:block, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ vsc(usc(environ), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=2)), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_6 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=1)), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    ### Identify the best model
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2, mod_henderson_3, mod_henderson_4, mod_henderson_5, mod_henderson_6), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2, mod_newtonrap_3, mod_newtonrap_4, mod_newtonrap_5, mod_newtonrap_6),
        verbose=verbose)
    list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=TRUE, verbose=verbose)
    if (is.na(list_u_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    if (verbose) {
        print(paste0("'Best fit': {", 
            "'algorithm': '", list_u_V_fitstats$algorithm, "', ", 
            "'model': {", 
                paste(paste0("'", names(list_u_V_fitstats$model), "': '", list_u_V_fitstats$model, "', "), collapse=" "),
            "}",
        "}"))
    }
    ### Extract BLUPs per environment
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod_newtonrap_5, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_6, mod_newtonrap=NA, verbose=verbose)
    df_BLUPs_GXE = fn_extract_gxe_breeding_values(list_u_V_fitstats)
    return(list(
        df_effects=df_BLUPs_GXE,
        loglik=list_u_V_fitstats$loglik,
        AIC=list_u_V_fitstats$AIC,
        BIC=list_u_V_fitstats$BIC,
        algorithm=list_u_V_fitstats$algorithm,
        model=list_u_V_fitstats$model))
}

#' Fit a GXE model to extract the best linear unbiased predictors for the genotype values,
#' assuming a spatially (SPAT) explicit design with row-by-column coordinates in a multiple environments.
#' We use this design when we expect significant heterogeneity in the each trial area, e.g. field trials with known nutrient gradients along one direction and moisture gradient along an orthogonal direction.
#' 
#' @param df data frame containing the model variables
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field in the data frame
#' @param env name of the environment field in the data frame
#' @param row name of the row field in the data frame
#' @param col name of the column field in the data frame
#' @param tolParInv diagonal loading of the variance-covariance matrix prior to inversion
#' @param skip_algo_based_on_dimensions if there are potentially more effects to be estimated than the number of observations then skip Henderson's fitting in favour of the Newton-Raphson transformations, and vice-versa
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values per environment (fields: id, env, BLUPs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the random effects
#' @examples
#' df = fn_simulate_gx1(design="spat")
#' out = fn_GXE_SPAT_BLUPs(df=df, trait="y", id="gen", env="env", row="row", col="col", 
#'      tolParInv=0.01, skip_algo_based_on_dimensions=TRUE, verbose=TRUE)
#' @export
fn_GXE_SPAT_BLUPs = function(df, trait="y", id="gen", env="env", row="row", col="col", tolParInv=0.01, skip_algo_based_on_dimensions=TRUE, verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gxe(design="spat")
    # trait="y"; id="gen"; env="env"; row="row"; col="col"; tolParInv=0.01; skip_algo_based_on_dimensions=TRUE; verbose=TRUE
    ### TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", env, "`"))))) {
        print(paste0("Error: env: ", env, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", row, "`"))))) {
        print(paste0("Error: row: ", row, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", col, "`"))))) {
        print(paste0("Error: col: ", col, " does not exist in the input dataframe."))
        return(NA)
    }
    ### Create consistently named response and explanatory variables for code readability below
    df_gxe_spat = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = eval(parse(text=paste0("df$", id))), 
        environ = eval(parse(text=paste0("df$", env))),
        row = eval(parse(text=paste0("df$", row))),
        col = eval(parse(text=paste0("df$", col)))
    )
    df_gxe_spat$row_factor = as.factor(df_gxe_spat$row)
    df_gxe_spat$col_factor = as.factor(df_gxe_spat$col)
    n_rows = nlevels(df_gxe_spat$row_factor)
    n_cols = nlevels(df_gxe_spat$col_factor)
    if (sum(!is.na(df_gxe_spat$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    mat_design_matrix_spat_gxe = stats::model.matrix(y ~ id:environ:row_factor:col_factor, data=df_gxe_spat)
    if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_spat_gxe) < ncol(mat_design_matrix_spat_gxe))) {
        ### if n < p (more coefficients to predict) then use Newton-Raphson transformations via mmer()
        skip_newtonrap = FALSE
        skip_henderson = TRUE
        mod_henderson_1 = mod_henderson_2 = mod_henderson_3 = mod_henderson_4 = mod_henderson_5 = mod_henderson_6 = NA
    } else if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_spat_gxe) >= ncol(mat_design_matrix_spat_gxe))) {
        ### if n >= p (less coefficients to predict) then use Henderson's mixed model equations via mmec()
        skip_newtonrap = TRUE
        skip_henderson = FALSE
        mod_newtonrap_1 = mod_newtonrap_2 = mod_newtonrap_3 = mod_newtonrap_4 = mod_newtonrap_5 = mod_newtonrap_6 = NA
    } else {
        skip_newtonrap = FALSE
        skip_henderson = FALSE
    }
    if (!skip_newtonrap) {
        if (verbose) {
            print("Fitting models via Newton-Raphson algorithm")
        }
        ###m Newton-Raphson models
        mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ environ:id + environ:row_factor:col_factor, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + environ:id + environ:row_factor:col_factor, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ vsr(usr(environ), id) + vsr(usr(environ), row_factor:col_factor),  rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id) + vsr(usr(environ), row_factor) + vsr(usr(environ), col_factor) + vsr(usr(environ), row_factor:col_factor),  rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ vsr(usr(environ), id) + sommer::spl2Da(x.coord=col, y.coord=row, nsegments=c(n_cols, n_rows), degree=c(3,3)),  rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_6 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id) + sommer::spl2Da(x.coord=col, y.coord=row, nsegments=c(n_cols, n_rows), degree=c(3,3)),  rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_7 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id) + vsr(usr(environ), row_factor) + vsr(usr(environ), col_factor) + sommer::spl2Da(x.coord=col, y.coord=row, nsegments=c(n_cols, n_rows), degree=c(3,3)),  rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    if (!skip_henderson) {
        if (verbose) {
            print("Fitting models via Henderson's mixed model equations")
        }
        ### Henderson's models
        mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ environ:id + environ:row_factor:col_factor, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + environ:id + environ:row_factor:col_factor, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ vsc(usc(environ), isc(id)) + vsc(usc(environ), isc(row_factor:col_factor)),  rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)) + vsc(usc(environ), isc(row_factor)) + vsc(usc(environ), isc(col_factor)) + vsc(usc(environ), isc(row_factor:col_factor)),  rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ vsc(usc(environ), isc(id)) + sommer::spl2Dc(x.coord=col, y.coord=row, nsegments=c(n_cols, n_rows), degree=c(3,3)),  rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_6 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)) + sommer::spl2Dc(x.coord=col, y.coord=row, nsegments=c(n_cols, n_rows), degree=c(3,3)),  rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_7 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)) + vsc(usc(environ), isc(row_factor)) + vsc(usc(environ), isc(col_factor)) + sommer::spl2Dc(x.coord=col, y.coord=row, nsegments=c(n_cols, n_rows), degree=c(3,3)),  rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_spat, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    ### Identify the best model
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2, mod_henderson_3, mod_henderson_4, mod_henderson_5, mod_henderson_6), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2, mod_newtonrap_3, mod_newtonrap_4, mod_newtonrap_5, mod_newtonrap_6),
        verbose=verbose)
    list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=TRUE, verbose=verbose)
    if (is.na(list_u_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    if (verbose) {
        print(paste0("'Best fit': {", 
            "'algorithm': '", list_u_V_fitstats$algorithm, "', ", 
            "'model': {", 
                paste(paste0("'", names(list_u_V_fitstats$model), "': '", list_u_V_fitstats$model, "', "), collapse=" "),
            "}",
        "}"))
    }
    ### Extract BLUPs per environment
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod_newtonrap_3, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_6, mod_newtonrap=NA, verbose=verbose)
    df_BLUPs_GXE = fn_extract_gxe_breeding_values(list_u_V_fitstats)
    return(list(
        df_effects=df_BLUPs_GXE,
        loglik=list_u_V_fitstats$loglik,
        AIC=list_u_V_fitstats$AIC,
        BIC=list_u_V_fitstats$BIC,
        algorithm=list_u_V_fitstats$algorithm,
        model=list_u_V_fitstats$model))
}

### Fit a 2-step GxE model