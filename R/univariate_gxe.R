### Assess genotype-by-environment (GXE) interactions by calculating pairwise trait correlations between environments.
### The trait values are simply aggregated using arithmetic mean and the pairwise correlations are calculated.
fn_assess_GXE = function(df, trait="y", id="gen", env="env", verbose=FALSE) {
    ### TEST #################################################################################
    # n = 100
    # G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
    # list_df_CORR = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=5, verbose=FALSE)
    # df=list_df_CORR$df; trait="y"; id="gen"; env="env"; verbose=TRUE
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

### Fit a GXE model to extract the best linear unbiased predictors for the genotype values,
### assuming a completely randomised design (CRD) in multiple environments.
### We use this design if per environment, we do not expect heterogeneity in the entire trial area, e.g. small trial and highly-controlled conditions.
fn_GXE_CRD_BLUPs = function(df, trait="y", id="gen", env="env", skip_algo_based_on_dimensions=TRUE, verbose=FALSE) {
    ### TEST #################################################################################
    # library(sommer)
    # n = 100
    # G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
    # list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=5, verbose=FALSE)
    # df=list_df_CORR_list_Y$df; trait="y"; id="gen"; env="env"; skip_algo_based_on_dimensions=TRUE; verbose=TRUE
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
        mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ environ:id, rcov= ~ units, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + environ:id, rcov= ~ units, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ vsr(usr(environ), id), rcov= ~ vsr(dsr(environ), units), data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id), rcov= ~ vsr(dsr(environ), units), data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=2)), id), rcov= ~ units, nIters=100, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_6 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=1)), id), rcov= ~ units, nIters=100, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    if (!skip_henderson) {
        if (verbose) {
            print("Fitting models via Henderson's mixed model equations")
        }
        ### Henderson's models
        mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ environ:id, rcov= ~ units, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + environ:id, rcov= ~ units, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ vsc(usc(environ), isc(id)), rcov= ~ vsc(dsc(environ), isc(units)), data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)), rcov= ~ vsc(dsc(environ), isc(units)), data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=2)), isc(id)), rcov= ~ units, nIters=100, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_6 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=1)), isc(id)), rcov= ~ units, nIters=100, data=df_gxe_crd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    ### Identify the best model
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2, mod_henderson_3, mod_henderson_4, mod_henderson_5, mod_henderson_6), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2, mod_newtonrap_3, mod_newtonrap_4, mod_newtonrap_5, mod_newtonrap_6),
        verbose=verbose)
    list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, verbose=verbose)
    if (is.na(list_u_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    if (verbose) {
        print(paste0("'Best fit': {", 
            "'algorithm': '", list_u_fitstats$algorithm, "', ", 
            "'model': {", 
                paste(paste0("'", names(list_u_fitstats$model), "': '", list_u_fitstats$model, "', "), collapse=" "),
            "}",
        "}"))
    }
    ### Extract BLUPs per environment
    formula_random_components = paste(as.character(list_u_fitstats$model$random), collapse="")
    vec_id = unique(df_gxe_crd$id)
    vec_environ = unique(df_gxe_crd$environ)
    vec_u_names = names(list_u_fitstats$u)


    # if ((list_u_fitstats$algorithm == "Newton-Raphson") & (formula_random_components == "~environ:id, rcov= ~ units")) {

    # } else if ((list_u_fitstats$algorithm == "Newton-Raphson") & (formula_random_components == "~id + environ:id, rcov= ~ units"))



    if (((sum(duplicated(vec_u_names)) > 0) & (sum(grepl("^PC", vec_u_names)) == 0)) | (length(vec_u_names)==length(vec_id))) {
        if (grepl("rrc\\(", formula_random_components)) {
            n_PC = length(vec_u_names) / length(vec_id)
            vec_u_names = paste0(rep(paste0("PC", 1:n_PC), each=length(vec_id)), ":id.y.", vec_u_names)
        } else if (length(vec_u_names)==length(vec_id)) {
            bool_idx_naked_entry_names = vec_u_names %in% vec_id
            vec_u_names[bool_idx_naked_entry_names] = paste0("id.y.id", vec_u_names[bool_idx_naked_entry_names])
            vec_u_names[!bool_idx_naked_entry_names] = paste0(rep(vec_environ, each=length(vec_id)), ":id.y.", vec_u_names[!bool_idx_naked_entry_names])
        } else {
            vec_u_names = paste0(rep(vec_environ, each=length(vec_id)), ":id.y.", vec_u_names)
        }
    }
    vec_u_names = gsub(paste0("environ:id.", trait, ".environ"), "", vec_u_names)
    if (!grepl("rrc\\(", formula_random_components)) {
        ### GxE variables per se (i.e. non-factor analytic)
        BLUPs = matrix(NA, nrow=length(vec_id), ncol=length(vec_environ))
        rownames(BLUPs) = vec_id
        colnames(BLUPs) = vec_environ
        for (i in 1:length(vec_id)) {
            for (j in 1:length(vec_environ)) {
                # i = 1; j = 2
                id_name = vec_id[i]
                environ_name = vec_environ[j]
                if (grepl("environ:id", formula_random_components)) {
                    pattern_gxe = paste0("^", environ_name)
                } else if (grepl("vsr\\(usr\\(environ\\), id\\)", formula_random_components) | 
                            grepl("vsr\\(dsr\\(environ\\), id\\)", formula_random_components) |
                            grepl("vsc\\(usc\\(environ\\), isc\\(id\\)\\)", formula_random_components) |
                            grepl("vsc\\(dsc\\(environ\\), isc\\(id\\)\\)", formula_random_components)) {
                    pattern_gxe = paste0("^", environ_name, ":id.y.")
                } else {
                    print("Error: add model specs here.")
                }
                idx_g = which(grepl(paste0("^id.y.id", id_name, "$"), vec_u_names)) 
                idx_gxe = which(grepl(id_name, vec_u_names) & grepl(pattern_gxe, vec_u_names))
                if (length(idx_g) == 0) {
                    BLUPs[i, j] = list_u_fitstats$u[idx_gxe]
                } else {
                    BLUPs[i, j] = list_u_fitstats$u[idx_g] + list_u_fitstats$u[idx_gxe]
                }
            }
        }
    } else {
        ### Factor analytic
        n_PC = as.numeric(unlist(strsplit(
            unlist(strsplit(gsub(" ", "", formula_random_components), "nPC="))[2], 
            ")"))[1])
        FA_gamma = with(df_gxe_crd, sommer::rrc(timevar=environ, idvar=id, response=y, nPC=n_PC, returnGamma = TRUE))$Gamma
        scores = matrix(NA, nrow=length(vec_id), ncol=n_PC)
        rownames(scores) = vec_id
        for (i in 1:length(vec_id)) {
            for (j in 1:n_PC) {
                # i = 1; j = 2
                id_name = vec_id[i]
                pc_name = paste0("PC", j)
                idx_g = which(grepl(paste0("^id.y.id", id_name, "$"), vec_u_names)) 
                idx_pc = which(grepl(id_name, vec_u_names) & grepl(paste0("^", pc_name, ":id.y."), vec_u_names))
                if (length(idx_g) == 0) {
                    scores[i, j] = list_u_fitstats$u[idx_pc]
                } else {
                    scores[i, j] = list_u_fitstats$u[idx_g] + list_u_fitstats$u[idx_pc]
                }
            }
        }
        BLUPs = scores %*% t(FA_gamma)
    }
    colnames(BLUPs) = paste0(trait, "_", colnames(BLUPs))
    df_BLUPs = eval(parse(text=paste0("data.frame(", id, "=rownames(BLUPs), BLUPs)"))) ### R converts dashes into dots
    rownames(df_BLUPs) = NULL
    colnames(df_BLUPs)[-1] = colnames(BLUPs) ### Revert to original column names of the BLUPs
    return(list(df_BLUPs=df_BLUPs, loglik=list_u_fitstats$loglik, AIC=list_u_fitstats$AIC, BIC=list_u_fitstats$BIC, algorithm=list_u_fitstats$algorithm, model=list_u_fitstats$model))
}

## Fit a GXE model to extract the best linear unbiased predictors for the genotype values,
## assuming a randomised complete block design (RCBD; each block is a complete replication of all genotypes) in multiple environments.
## We use this design when we expect heterogeneity along a single direction in the trial area, e.g. field trial along a slope.
fn_GXE_RCBD_BLUPs = function(df, trait="y", id="gen", env="env", block="rep", skip_algo_based_on_dimensions=TRUE, verbose=FALSE) {
    ### TEST #################################################################################
    # library(sommer)
    # n = 100
    # n_blocks = 5
    # G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
    # list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=n_blocks, verbose=FALSE)
    # df=list_df_CORR_list_Y$df; trait="y"; id="gen"; env="env"; skip_algo_based_on_dimensions=TRUE; verbose=TRUE
    # ### Simulate block effects
    # for (j in 1:n_blocks) {
    #     idx = which(df$rep == as.character(j))
    #     df$y[idx] = df$y[idx] + stats::rnorm(1)
    # }
    # trait="y"; id="gen"; env="env"; block="rep"; skip_algo_based_on_dimensions=TRUE; verbose=TRUE
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
    df_gxe_rcbd = data.frame(
        y = eval(parse(text=paste0("df$", trait))), 
        id = eval(parse(text=paste0("df$", id))), 
        environ = eval(parse(text=paste0("df$", env))),
        block = eval(parse(text=paste0("df$", block)))
    )
    if (sum(!is.na(df_gxe_rcbd$y)) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    mat_design_matrix_rcbd_gxe = stats::model.matrix(y ~ id:environ:block, data=df_gxe_rcbd)
    if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_rcbd_gxe) < ncol(mat_design_matrix_rcbd_gxe))) {
        ### if n < p (more coefficients to predict) then use Newton-Raphson transformations via mmer()
        skip_newtonrap = FALSE
        skip_henderson = TRUE
        mod_henderson_1 = mod_henderson_2 = mod_henderson_3 = mod_henderson_4 = mod_henderson_5 = mod_henderson_6 = NA
    } else if (skip_algo_based_on_dimensions & (nrow(mat_design_matrix_rcbd_gxe) >= ncol(mat_design_matrix_rcbd_gxe))) {
        ### if n >= p (less coefficients to predict) then use Henderson's mixed model equations via mmec()
        skip_newtonrap = TRUE
        skip_henderson = FALSE
        mod_newtonrap_1 = mod_newtonrap_2 = mod_newtonrap_3 = mod_newtonrap_4 = mod_newtonrap_5 = mod_newtonrap_6 = NA
    } else {
        skip_newtonrap = FALSE
        skip_henderson = FALSE
    }
    tolParInv=0.01
    if (!skip_newtonrap) {
        if (verbose) {
            print("Fitting models via Newton-Raphson algorithm")
        }
        ###m Newton-Raphson models
        mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ environ:id + environ:block, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + environ:id + environ:block, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ vsr(usr(environ), id) + vsr(usr(environ), block), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id) + vsr(usr(environ), block), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=2)), id) + vsr(usr(environ), block), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_6 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=1)), id) + vsr(usr(environ), block), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    if (!skip_henderson) {
        if (verbose) {
            print("Fitting models via Henderson's mixed model equations")
        }
        ### Henderson's models
        mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ environ:id + environ:block, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + environ:id, rcov= ~ units, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ vsc(usc(environ), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=2)), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_6 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=1)), isc(id)) + vsc(usc(environ), isc(block)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df_gxe_rcbd, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
    }
    ### Identify the best model
    list_best_mods = fn_find_best_fit_within_algo(
        list_mod_henderson=list(mod_henderson_1, mod_henderson_2, mod_henderson_3, mod_henderson_4, mod_henderson_5, mod_henderson_6), 
        list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2, mod_newtonrap_3, mod_newtonrap_4, mod_newtonrap_5, mod_newtonrap_6),
        verbose=verbose)
    list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=list_best_mods$mod_henderson, mod_newtonrap=list_best_mods$mod_newtonrap, verbose=verbose)
    if (is.na(list_u_fitstats[1])) {
        print("Error: failed to fit linear mixed models.")
        return(NA)
    }
    if (verbose) {
        print(paste0("'Best fit': {", 
            "'algorithm': '", list_u_fitstats$algorithm, "', ", 
            "'model': {", 
                paste(paste0("'", names(list_u_fitstats$model), "': '", list_u_fitstats$model, "', "), collapse=" "),
            "}",
        "}"))
    }
    ### Extract BLUPs per environment
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_newtonrap=mod_newtonrap_1, mod_henderson=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_newtonrap=mod_newtonrap_2, mod_henderson=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_newtonrap=mod_newtonrap_3, mod_henderson=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_newtonrap=mod_newtonrap_4, mod_henderson=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_newtonrap=mod_newtonrap_5, mod_henderson=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_newtonrap=mod_newtonrap_6, mod_henderson=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_1, mod_newtonrap=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_2, mod_newtonrap=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_3, mod_newtonrap=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_4, mod_newtonrap=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_5, mod_newtonrap=NA)
    # list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson_6, mod_newtonrap=NA)



    ### !!!!NEED A BETTER NAMES PARSING!!!!



    formula_random_components = paste(as.character(list_u_fitstats$model$random), collapse="")
    vec_id = unique(df_gxe_rcbd$id)
    vec_environ = unique(df_gxe_rcbd$environ)
    vec_u_names = names(list_u_fitstats$u)
    # vec_u_names[!grepl("block", vec_u_names)] = NA ### Remove all random effects with block effects as we are only interested in the id and environmental effects
    if (((sum(duplicated(vec_u_names)) > 0) & (sum(grepl("^PC", vec_u_names)) == 0)) | (length(vec_u_names)==length(vec_id))) {
        if (grepl("rrc\\(", formula_random_components)) {
            n_PC = length(vec_u_names) / length(vec_id)
            vec_u_names = paste0(rep(paste0("PC", 1:n_PC), each=length(vec_id)), ":id.y.", vec_u_names)
        } else if (length(vec_u_names)==length(vec_id)) {
            bool_idx_naked_entry_names = vec_u_names %in% vec_id
            vec_u_names[bool_idx_naked_entry_names] = paste0("id.y.id", vec_u_names[bool_idx_naked_entry_names])
            vec_u_names[!bool_idx_naked_entry_names] = paste0(rep(vec_environ, each=length(vec_id)), ":id.y.", vec_u_names[!bool_idx_naked_entry_names])
        } else {
            vec_u_names = paste0(rep(vec_environ, each=length(vec_id)), ":id.y.", vec_u_names)
        }
    }
    vec_u_names = gsub(paste0("environ:id.", trait, ".environ"), "", vec_u_names)
    if (!grepl("rrc\\(", formula_random_components)) {
        ### GxE variables per se (i.e. non-factor analytic)
        BLUPs = matrix(NA, nrow=length(vec_id), ncol=length(vec_environ))
        rownames(BLUPs) = vec_id
        colnames(BLUPs) = vec_environ
        for (i in 1:length(vec_id)) {
            for (j in 1:length(vec_environ)) {
                # i = 1; j = 2
                id_name = vec_id[i]
                environ_name = vec_environ[j]
                if (grepl("environ:id", formula_random_components)) {
                    pattern_gxe = paste0("^", environ_name)
                } else if (grepl("vsr\\(usr\\(environ\\), id\\)", formula_random_components) | 
                            grepl("vsr\\(dsr\\(environ\\), id\\)", formula_random_components) |
                            grepl("vsc\\(usc\\(environ\\), isc\\(id\\)\\)", formula_random_components) |
                            grepl("vsc\\(dsc\\(environ\\), isc\\(id\\)\\)", formula_random_components)) {
                    pattern_gxe = paste0("^", environ_name, ":id.y.")
                } else {
                    print("Error: add model specs here.")
                }
                idx_g = which(grepl(paste0("^id.y.id", id_name, "$"), vec_u_names)) 
                idx_gxe = which(grepl(id_name, vec_u_names) & grepl(pattern_gxe, vec_u_names))
                if (length(idx_g) == 0) {
                    BLUPs[i, j] = list_u_fitstats$u[idx_gxe]
                } else {
                    BLUPs[i, j] = list_u_fitstats$u[idx_g] + list_u_fitstats$u[idx_gxe]
                }
            }
        }
    } else {
        ### Factor analytic
        n_PC = as.numeric(unlist(strsplit(
            unlist(strsplit(gsub(" ", "", formula_random_components), "nPC="))[2], 
            ")"))[1])
        FA_gamma = with(df_gxe_rcbd, sommer::rrc(timevar=environ, idvar=id, response=y, nPC=n_PC, returnGamma = TRUE))$Gamma
        scores = matrix(NA, nrow=length(vec_id), ncol=n_PC)
        rownames(scores) = vec_id
        for (i in 1:length(vec_id)) {
            for (j in 1:n_PC) {
                # i = 1; j = 2
                id_name = vec_id[i]
                pc_name = paste0("PC", j)
                idx_g = which(grepl(paste0("^id.y.id", id_name, "$"), vec_u_names)) 
                idx_pc = which(grepl(id_name, vec_u_names) & grepl(paste0("^", pc_name, ":id.y."), vec_u_names))
                if (length(idx_g) == 0) {
                    scores[i, j] = list_u_fitstats$u[idx_pc]
                } else {
                    scores[i, j] = list_u_fitstats$u[idx_g] + list_u_fitstats$u[idx_pc]
                }
            }
        }
        BLUPs = scores %*% t(FA_gamma)
    }
    colnames(BLUPs) = paste0(trait, "_", colnames(BLUPs))
    df_BLUPs = eval(parse(text=paste0("data.frame(", id, "=rownames(BLUPs), BLUPs)"))) ### R converts dashes into dots
    rownames(df_BLUPs) = NULL
    colnames(df_BLUPs)[-1] = colnames(BLUPs) ### Revert to original column names of the BLUPs
    return(list(df_BLUPs=df_BLUPs, loglik=list_u_fitstats$loglik, AIC=list_u_fitstats$AIC, BIC=list_u_fitstats$BIC, algorithm=list_u_fitstats$algorithm, model=list_u_fitstats$model))
}
