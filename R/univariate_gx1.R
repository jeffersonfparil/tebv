#' Fit mixed models to extract the best linear unbiased predictors for the genotype values,
#' assuming a completely randomised design (CRD) in a single environment.
#' We use this design when we do not expect heterogeneity in the entire trial area, e.g. small trial and highly-controlled conditions.
#' 
#' @param df data frame containing the model variables
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field (main explanatory variable) in the data frame
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values (fields: id, BLUPs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the random effects
#' df_fixed_effects: data frame of fixed factor names, and their corresponding effects
#' @examples
#' df = fn_simulate_gx1(design="crd")
#' out = fn_GX1_CRD_BLUPs(df=df, trait="y", id="gen", verbose=TRUE)
#' @export
fn_GX1_CRD_BLUPs = function(df, trait="y", id="gen", verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gx1(design="crd")
    # trait="y"; id="gen"; verbose=TRUE
    ## TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    ### Create consistently named response and explanatory variables for code readability below
    df_crd = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = as.factor(eval(parse(text=paste0("df$", id))))
    )
    if (sum(!is.na(df_crd$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    ### Fit
    mod_henderson = tryCatch(sommer::mmec(y ~ 1, random = ~ id, rcov= ~ units, data=df_crd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap = tryCatch(sommer::mmer(y ~ 1, random = ~ id, rcov= ~ units, data=df_crd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson, mod_newtonrap=mod_newtonrap, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson, mod_newtonrap=NA, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod_newtonrap, verbose=verbose)
    if (is.na(list_u_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    u = list_u_V_fitstats$u
    df_BLUPs = eval(parse(text=paste0("data.frame(", id, "=names(u), ", trait, "=u)")))
    rownames(df_BLUPs) = NULL
    return(list(
        df_effects=df_BLUPs, 
        loglik=list_u_V_fitstats$loglik,
        AIC=list_u_V_fitstats$AIC,
        BIC=list_u_V_fitstats$BIC,
        algorithm=list_u_V_fitstats$algorithm,
        model=list_u_V_fitstats$model,
        V=list_u_V_fitstats$V))
}

#' Fit linear models to extract the best linear unbiased estimators for the genotype values,
#' assuming a completely randomised design (CRD) in a single environment.
#' We use this design when we do not expect heterogeneity in the entire trial area, e.g. small trial and highly-controlled conditions.
#' 
#' @param df data frame containing the model variables
#' @param control_id name of the control genotype you wish to set as the the intercept from which all other entries where compared against
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field (main explanatory variable) in the data frame
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values (fields: id, BLUEs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the fixed effects
#' df_fixed_effects: data frame of fixed factor names, and their corresponding effects
#' @examples
#' df = fn_simulate_gx1(design="crd")
#' out = fn_GX1_CRD_BLUEs(df=df, control_id=unique(df$gen)[2], trait="y", id="gen", verbose=TRUE)
#' @export
fn_GX1_CRD_BLUEs = function(df, control_id, trait="y", id="gen", verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gx1(design="crd")
    # control_id=unique(df$gen)[2]
    # trait="y"; id="gen"; verbose=TRUE
    ## TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    ### Create consistently named response and explanatory variables for code readability below
    df_crd = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = as.factor(eval(parse(text=paste0("df$", id))))
    )
    if (sum(!is.na(df_crd$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    ### Set control as the intercept so that all the BLUEs are relative to the control genotype
    df_crd$id = stats::relevel(df_crd$id, ref=control_id)
    ### Fit
    mod_henderson = tryCatch(sommer::mmec(y ~ 1 + id, rcov= ~ units, data=df_crd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap = tryCatch(sommer::mmer(y ~ 1 + id, rcov= ~ units, data=df_crd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    ### Extract BLUEs + intercept effect including the intercept renamed as the control genotype
    list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson, mod_newtonrap=mod_newtonrap, extract_BLUPs=FALSE, verbose=verbose)
    # list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson, mod_newtonrap=NA, extract_BLUPs=FALSE, verbose=verbose)
    # list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod_newtonrap, extract_BLUPs=FALSE, verbose=verbose)
    if (is.na(list_b_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    b = list_b_V_fitstats$b
    df_BLUEs = eval(parse(text=paste0("data.frame(", id, "=names(b), ", trait, "=b)")))
    rownames(df_BLUEs) = NULL
    return(list(
        df_effects=df_BLUEs, 
        loglik=list_b_V_fitstats$loglik,
        AIC=list_b_V_fitstats$AIC,
        BIC=list_b_V_fitstats$BIC,
        algorithm=list_b_V_fitstats$algorithm,
        model=list_b_V_fitstats$model,
        V=list_b_V_fitstats$V))
}

#' Fit mixed models to extract the best linear unbiased predictors for the genotype values,
#' assuming a randomised complete/incomplete block design (RBD; each block is a complete replication of all genotypes) in a single environment.
#' We use this design when we expect heterogeneity along a single direction in the trial area, e.g. field trial along a slope.
#' 
#' @param df data frame containing the model variables
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field (main explanatory variable) in the data frame
#' @param block name of the complete or incomplete block field in the data frame
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values (fields: id, BLUPs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the random effects
#' df_fixed_effects: data frame of fixed factor names, and their corresponding effects
#' @examples
#' df = fn_simulate_gx1(design="rbd")
#' out = fn_GX1_RBD_BLUPs(df=df, trait="y", id="gen", block="rep", verbose=TRUE)
#' @export
fn_GX1_RBD_BLUPs = function(df, trait="y", id="gen", block="rep", verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gx1(design="rbd")
    # trait="y"; id="gen"; block="rep"; verbose=TRUE
    ## TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", block, "`"))))) {
        print(paste0("Error: block: ", block, " does not exist in the input dataframe."))
        return(NA)
    }
    ### Create consistently named response and explanatory variables for code readability below
    df_rbd = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = as.factor(eval(parse(text=paste0("df$", id)))), 
        block = as.factor(eval(parse(text=paste0("df$", block))))
    )
    if (sum(!is.na(df_rbd$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + block, random = ~ id, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1, random = ~ id + block, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + block, random = ~ id, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1, random = ~ id + block, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    ### Find the best model per algorithm
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2),
        verbose=verbose)
    list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=NA, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=list_best_mods$mod_newtonrap, verbose=verbose)
    if (is.na(list_u_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    idx = which(names(list_u_V_fitstats$u) %in% unique(df_rbd$id))
    u = list_u_V_fitstats$u[idx]
    df_BLUPs = eval(parse(text=paste0("data.frame(", id, "=names(u), ", trait, "=u)")))
    rownames(df_BLUPs) = NULL
    return(list(
        df_effects=df_BLUPs, 
        loglik=list_u_V_fitstats$loglik,
        AIC=list_u_V_fitstats$AIC,
        BIC=list_u_V_fitstats$BIC,
        algorithm=list_u_V_fitstats$algorithm,
        model=list_u_V_fitstats$model,
        V=list_u_V_fitstats$V))
}

#' Fit mixed models to extract the best linear unbiased estimators for the genotype values,
#' assuming a randomised complete/incomplete block design (RBD; each block is a complete replication of all genotypes) in a single environment.
#' We use this design when we expect heterogeneity along a single direction in the trial area, e.g. field trial along a slope.
#' 
#' @param df data frame containing the model variables
#' @param control_id name of the control genotype you wish to set as the the intercept from which all other entries where compared against
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field (main explanatory variable) in the data frame
#' @param block name of the complete or incomplete block field in the data frame
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values (fields: id, BLUEs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the fixed effects
#' df_fixed_effects: data frame of fixed factor names, and their corresponding effects
#' @examples
#' df = fn_simulate_gx1(design="rbd")
#' out = fn_GX1_RBD_BLUEs(df=df, control_id=unique(df$gen)[3], trait="y", id="gen", 
#'      block="rep", verbose=TRUE)
#' @export
fn_GX1_RBD_BLUEs = function(df, control_id, trait="y", id="gen", block="rep", verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gx1(design="rbd")
    # control_id=unique(df$gen)[2]
    # trait="y"; id="gen"; block="rep"; verbose=TRUE
    ## TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", block, "`"))))) {
        print(paste0("Error: block: ", block, " does not exist in the input dataframe."))
        return(NA)
    }
    ### Create consistently named response and explanatory variables for code readability below
    df_rbd = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = as.factor(eval(parse(text=paste0("df$", id)))), 
        block = as.factor(eval(parse(text=paste0("df$", block))))
    )
    if (sum(!is.na(df_rbd$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    ### Set control as the intercept so that all the BLUEs are relative to the control genotype
    df_rbd$id = stats::relevel(df_rbd$id, ref=control_id)
    ### Fit
    mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + id, random = ~ block, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + id + block, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + id, random = ~ block, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + id + block, rcov= ~ units, data=df_rbd, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    ### Find the best model per algorithm
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2),
        verbose=verbose)
    ### Extract BLUEs + intercept effect including the intercept renamed as the control genotype
    list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=FALSE, verbose=verbose)
    # list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=NA, extract_BLUPs=FALSE, verbose=verbose)
    # list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=FALSE, verbose=verbose)
    if (is.na(list_b_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    idx = which(names(list_b_V_fitstats$b) %in% unique(df_rbd$id))
    b = list_b_V_fitstats$b[idx]
    df_BLUEs = eval(parse(text=paste0("data.frame(", id, "=names(b), ", trait, "=b)")))
    rownames(df_BLUEs) = NULL
    return(list(
        df_effects=df_BLUEs, 
        loglik=list_b_V_fitstats$loglik,
        AIC=list_b_V_fitstats$AIC,
        BIC=list_b_V_fitstats$BIC,
        algorithm=list_b_V_fitstats$algorithm,
        model=list_b_V_fitstats$model,
        V=list_b_V_fitstats$V))
}

#' Fit mixed models to extract the best linear unbiased predictors for the genotype values,
#' assuming a spatially (SPAT) explicit design with row-by-column coordinates in a single environment.
#' We use this design when we expect significant heterogeneity in the entire trial area, e.g. field trials with known nutrient gradients along one direction and moisture gradient along an orthogonal direction.
#' 
#' @param df data frame containing the model variables
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field (main explanatory variable) in the data frame
#' @param row name of the row field in the data frame
#' @param col name of the column field in the data frame
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values (fields: id, BLUPs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the random effects
#' df_fixed_effects: data frame of fixed factor names, and their corresponding effects
#' @examples
#' df = fn_simulate_gx1(design="spat")
#' out = fn_GX1_SPAT_BLUPs(df=df, trait="y", id="gen", row="row", col="col", verbose=TRUE)
#' @export
fn_GX1_SPAT_BLUPs = function(df, trait="y", id="gen", row="row", col="col", verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gx1(design="spat")
    # trait="y"; id="gen"; row="row"; col="col"; verbose=TRUE
    ## TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
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
    df_spat = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = as.factor(eval(parse(text=paste0("df$", id)))), 
        row = eval(parse(text=paste0("df$", row))), 
        col = eval(parse(text=paste0("df$", col)))
    )
    if (sum(!is.na(df_spat$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    if (!is.numeric(df_spat$row)) {
        print(paste0("Error: row: ", row, " is not numeric. Numeric row is required for spatial modelling."))
        return(NA)
    }
    if (!is.numeric(df_spat$col)) {
        print(paste0("Error: col: ", col, " is not numeric. Numeric column is required for spatial modelling."))
        return(NA)
    }
    df_spat$row_factor = as.factor(df_spat$row)
    df_spat$col_factor = as.factor(df_spat$col)
    n_rows = nlevels(df_spat$row_factor)
    n_cols = nlevels(df_spat$col_factor)
    ### Fit 5 models per algorithm
    mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1, random = ~id + row + col + sommer::spl2Dc(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1, random = ~id + row_factor + col_factor + sommer::spl2Dc(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + row_factor + col_factor, random = ~id + sommer::spl2Dc(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + row + col + row:col, random = ~id + sommer::spl2Dc(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + row + col + row:col, random = ~id + row_factor + col_factor + sommer::spl2Dc(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1, random = ~id + row + col + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1, random = ~id + row_factor + col_factor + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + row_factor + col_factor, random = ~id + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + row + col + row:col, random = ~id + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + row + col + row:col, random = ~id + row_factor + col_factor + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    ### Find the best model per algorithm
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2, mod_henderson_3, mod_henderson_4, mod_henderson_5), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2, mod_newtonrap_3, mod_newtonrap_4, mod_newtonrap_5),
        verbose=verbose)
    list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=TRUE, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=NA, verbose=verbose)
    # list_u_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=list_best_mods$mod_newtonrap, verbose=verbose)
    if (is.na(list_u_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    idx = which(names(list_u_V_fitstats$u) %in% unique(df_spat$id))
    u = list_u_V_fitstats$u[idx]
    df_BLUPs = eval(parse(text=paste0("data.frame(", id, "=names(u), ", trait, "=u)")))
    rownames(df_BLUPs) = NULL
    return(list(
        df_effects=df_BLUPs, 
        loglik=list_u_V_fitstats$loglik,
        AIC=list_u_V_fitstats$AIC,
        BIC=list_u_V_fitstats$BIC,
        algorithm=list_u_V_fitstats$algorithm,
        model=list_u_V_fitstats$model,
        V=list_u_V_fitstats$V))
}

#' Fit mixed models to extract the best linear unbiased estimators for the genotype values,
#' assuming a spatially (SPAT) explicit design with row-by-column coordinates in a single environment.
#' We use this design when we expect significant heterogeneity in the entire trial area, e.g. field trials with known nutrient gradients along one direction and moisture gradient along an orthogonal direction.
#' 
#' @param df data frame containing the model variables
#' @param control_id name of the control genotype you wish to set as the the intercept from which all other entries where compared against
#' @param trait name of the continuous numeric response variable in the data frame
#' @param id name of the entry field (main explanatory variable) in the data frame
#' @param row name of the row field in the data frame
#' @param col name of the column field in the data frame
#' @param verbose show model fitting messages?
#' @returns
#' df_effects: data frame trial-estimated breeding values (fields: id, BLUEs)
#' loglik: log-likelihood of the best-fitting model
#' AIC:  Akaike information criterion (prediction error estimator) of the best-fitting model
#' BIC: Bayesian information criterion (another prediction error estimator) of the best-fitting model
#' algorithm: model fitting algorithm used, i.e. Henderson's (mmec) or Newton-Raphson-transformations (mmer)
#' model: Specification of the best-fitting linear model
#' V: Variance-covariance component of the fixed effects
#' df_fixed_effects: data frame of fixed factor names, and their corresponding effects
#' @examples
#' df = fn_simulate_gx1(design="spat")
#' out = fn_GX1_SPAT_BLUEs(df=df, control_id=unique(df$gen)[4], trait="y", id="gen", 
#'      row="row", col="col", verbose=TRUE)
#' @export
fn_GX1_SPAT_BLUEs = function(df, control_id, trait="y", id="gen", row="row", col="col", verbose=FALSE) {
    ### TEST #################################################################################
    # source("helpers.R")
    # library(sommer)
    # df = fn_simulate_gx1(design="spat")
    # control_id=unique(df$gen)[2]
    # trait="y"; id="gen"; row="row"; col="col"; verbose=TRUE
    ## TEST #################################################################################
    if (is.null(eval(parse(text=paste0("df$`", trait, "`"))))) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    if (is.null(eval(parse(text=paste0("df$`", id, "`"))))) {
        print(paste0("Error: id: ", id, " does not exist in the input dataframe."))
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
    df_spat = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = as.factor(eval(parse(text=paste0("df$", id)))), 
        row = eval(parse(text=paste0("df$", row))), 
        col = eval(parse(text=paste0("df$", col)))
    )
    if (sum(!is.na(df_spat$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    if (!is.numeric(df_spat$row)) {
        print(paste0("Error: row: ", row, " is not numeric. Numeric row is required for spatial modelling."))
        return(NA)
    }
    if (!is.numeric(df_spat$col)) {
        print(paste0("Error: col: ", col, " is not numeric. Numeric column is required for spatial modelling."))
        return(NA)
    }
    ### Set control as the intercept so that all the BLUEs are relative to the control genotype
    df_spat$id = stats::relevel(df_spat$id, ref=control_id)
    ### Fit
    df_spat$row_factor = as.factor(df_spat$row)
    df_spat$col_factor = as.factor(df_spat$col)
    mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + id, random = ~row + col + row_factor:col_factor, rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + id, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + id + row + col, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + id + row + col + row:col, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + id + row_factor + col_factor + row_factor:col_factor, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    n_rows = nlevels(df_spat$row_factor)
    n_cols = nlevels(df_spat$col_factor)
    mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + id, random = ~row + col + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + id, random = ~row_factor + col_factor + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + id + row_factor + col_factor, random = ~sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + id + row + col + row:col, random = ~sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + id + row + col + row:col, random = ~row_factor + col_factor + sommer::spl2Da(x.coord=df_spat$col, y.coord=df_spat$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df_spat, dateWarning=FALSE, verbose=verbose),
        error=function(e){NA})
    ### Find the best model per algorithm
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2, mod_henderson_3, mod_henderson_4, mod_henderson_5), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2, mod_newtonrap_3, mod_newtonrap_4, mod_newtonrap_5),
        verbose=verbose)
    list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=FALSE, verbose=verbose)
    # list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=NA, extract_BLUPs=FALSE, verbose=verbose)
    # list_b_V_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=list_best_mods$mod_newtonrap, extract_BLUPs=FALSE, verbose=verbose)
    if (is.na(list_b_V_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    idx = which(names(list_b_V_fitstats$b) %in% unique(df_spat$id))
    b = list_b_V_fitstats$b[idx]
    df_BLUEs = eval(parse(text=paste0("data.frame(", id, "=names(b), ", trait, "=b)")))
    rownames(df_BLUEs) = NULL
    return(list(
        df_effects=df_BLUEs, 
        loglik=list_b_V_fitstats$loglik,
        AIC=list_b_V_fitstats$AIC,
        BIC=list_b_V_fitstats$BIC,
        algorithm=list_b_V_fitstats$algorithm,
        model=list_b_V_fitstats$model,
        V=list_b_V_fitstats$V))
}
