### We have decided to divide TOTAL_CBD such that we have 2 additional traits - 3 traits in total:
###     [1] original trait values, 
###     [2] lower range where values above the 1 becomes NA, and
###     [3] upper range where values below the 1 becomes NA.
fn_divide_into_lower_and_upper = function(df, trait="y", threshold=10, verbose=FALSE) {
    ### TEST #################################################################################
    # G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=2, verbose=TRUE)
    # list_df_CORR = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=5, verbose=TRUE)
    # txtplot::txtdensity(list_df_CORR$df$y)
    # ### Create a bimodal phenotype distribution by adding 5 to the phenotype values of 20% of the genotypes
    # vec_gen_add_10 = sort(sample(unique(list_df_CORR$df$gen), size=20, replace=FALSE))
    # idx = which(list_df_CORR$df$gen %in% vec_gen_add_10)
    # list_df_CORR$df$y[idx] = list_df_CORR$df$y[idx] + 10.0
    # txtplot::txtdensity(list_df_CORR$df$y)
    # trait="y"; threshold=10; verbose=TRUE
    ### TEST #################################################################################
    y = eval(parse(text=paste0("df$`", trait, "`")))
    if (is.null(y)) {
        print(paste0("Error: trait: ", trait, " does not exist in the input dataframe."))
        return(NA)
    }
    y = y[!is.na(y)]
    if (length(y) == 0) {
        print(paste0("Error: all data are missing for trait: ", trait, "."))
        return(NA)
    }
    idx_lower = which(y<threshold)
    idx_upper = which(y>=threshold)
    if (length(idx_lower) == 0) {
        print(paste0("Warning: there is no data below the threshold: ", threshold, " for the trait: ", trait, "."))
    }
    if (length(idx_upper) == 0) {
        print(paste0("Warning: there is no data above/equal the threshold: ", threshold, " for the trait: ", trait, "."))
    }
    ### Set value beyond the threshold to missing
    y_lower = y; y_lower[idx_upper] = NA
    y_upper = y; y_upper[idx_lower] = NA
    if (verbose) {
        print("################################")
        print(paste0("All (n=", sum(!is.na(y)), "):"))
        txtplot::txtdensity(y[!is.na(y)])
        print(paste0("Lower (n=", sum(!is.na(y_lower)), "):"))
        txtplot::txtdensity(y_lower[!is.na(y_lower)])
        print(paste0("Upper (n=", sum(!is.na(y_upper)), "):"))
        txtplot::txtdensity(y_upper[!is.na(y_upper)])
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    }
    ### Output the dataframe
    eval(parse(text=paste0("df$`", trait, "_LOWER` = y_lower")))
    eval(parse(text=paste0("df$`", trait, "_UPPER` = y_upper")))
    return(df)
}

### Find the best Henderson's and Newton-Raphson models when there are multiple models per algorithm.
### The model with higher log-likelihood, lower AIC, and lower BIC on average will be selected.
fn_find_best_fit_within_algo = function(list_mod_henderson, list_mod_newtonrap, verbose=FALSE) {
    ### TEST #################################################################################
    # n_alleles = 2
    # pheno_reps = 3
    # G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=n_alleles, verbose=FALSE)
    # list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, h2=0.75, pheno_reps=pheno_reps, verbose=FALSE)
    # Y = list_Y_b_E_b_epi$Y
    # df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=pheno_reps))
    # ### True genotype values
    # df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
    # ### Fit using mmec and mmer
    # mod_henderson_1 = tryCatch(suppressWarnings(sommer::mmec(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
    # mod_henderson_2 = tryCatch(suppressWarnings(sommer::mmec(y ~ 0, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
    # mod_newtonrap_1 = tryCatch(suppressWarnings(sommer::mmer(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
    # mod_newtonrap_2 = tryCatch(suppressWarnings(sommer::mmer(y ~ 0, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
    # list_mod_henderson=list(mod_henderson_1, mod_henderson_2)
    # list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2)
    # verbose = TRUE
    ### TEST #################################################################################
    if (!is.null(list_mod_henderson) & !is.null(list_mod_newtonrap)) {
        if (length(list_mod_henderson) != length(list_mod_newtonrap)) {
            print("Error: the number of models fitted using Henderson's mixed model equations and Newton-Raphson transformations do not match.")
            return(NA)
        }
        n = length(list_mod_henderson)
        vec_algo = c("henderson", "newtonrap")
    } else {
        if (!is.null(list_mod_henderson)) {
            n = length(list_mod_henderson)
            vec_algo = c("henderson")
            mod_newtonrap = NA
        } else {
            n = length(list_mod_newtonrap)
            vec_algo = c("newtonrap")
            mod_henderson = NA
        }
    }
    for (algo in vec_algo) {
        # algo = "henderson"
        # algo = "newtonrap"
        best_fitstats = NULL
        model_number_selected = 1
        for (i in 1:n) {
            # i = 1
            mod = eval(parse(text=paste0("list_mod_", algo, "[[", i, "]]")))
            if (!is.na(mod[1])) {
                tmp_fitstats = suppressWarnings(summary(mod)$logo)
                if (verbose) {
                    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                    if (algo == "henderson") {
                        print(paste0("ALGO: ", algo, "; fixed: ", paste(mod$args$fixed, collapse=" "), "; random: ", paste(mod$args$random, collapse=" ")))
                    } else {
                        print(paste0("ALGO: ", algo, "; fixed: ", paste(mod$call$fixed, collapse=" "), "; random: ", paste(mod$call$random, collapse=" ")))
                    }
                    print(tmp_fitstats)
                }
                if (is.null(best_fitstats)) {
                    best_fitstats = tmp_fitstats
                    model_number_selected = i
                } else {
                    fitstats_comparison = mean(c(
                        (tmp_fitstats$logLik > best_fitstats$logLik), 
                        (tmp_fitstats$AIC < best_fitstats$AIC), 
                        (tmp_fitstats$BIC < best_fitstats$BIC)
                    ))
                    if (fitstats_comparison >= 0.5) {
                        best_fitstats = tmp_fitstats
                        model_number_selected = i
                    }
                }
            }
        }
        eval(parse(text=paste0("mod_", algo, " = list_mod_", algo, "[[", model_number_selected, "]]")))
    }
    return(list(mod_henderson=mod_henderson, mod_newtonrap=mod_newtonrap))
}

### Extracting BLUPs from the Henderson's mixed model equations-based fit or Newton-Raphson's transformations-based  (mmec vs mmer)
### depending on which has higher log-likelihood, lower AIC, and lower BIC on average.
fn_henderson_vs_newtonraphson_fit = function(mod_henderson, mod_newtonrap, verbose=FALSE) {
    ### TEST #################################################################################
    # G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=2, verbose=TRUE)
    # list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=2, dist_effects="norm", n_effects=100, h2=0.75, pheno_reps=3, verbose=TRUE)
    # Y = list_Y_b_E_b_epi$Y
    # df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=3))
    # mod_henderson = tryCatch(sommer::mmec(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=verbose), error=function(e){NA})
    # mod_newtonrap = tryCatch(sommer::mmer(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=verbose), error=function(e){NA})
    ### TEST #################################################################################
    if (is.na(mod_henderson[1]) & is.na(mod_newtonrap[1])) {
        return(NA)
    }
    ### Select the model that is non-missing or have better fit
    if (is.na(mod_henderson[1])) {
        fitstats_comparison = 0.0
        fitstats_newtonrap = summary(mod_newtonrap)$logo
    } else if (is.na(mod_newtonrap[1])) {
        fitstats_comparison = 1.0
        fitstats_henderson = summary(mod_henderson)$logo
    } else {
        fitstats_henderson = summary(mod_henderson)$logo
        fitstats_newtonrap = summary(mod_newtonrap)$logo
        if (verbose) {
            print("Henderson fit:")
            print(fitstats_henderson)
            print("Newton-Raphson fit:")
            print(fitstats_newtonrap)
        }
        fitstats_comparison = mean(c(
            (fitstats_henderson$logLik > fitstats_newtonrap$logLik), 
            (fitstats_henderson$AIC < fitstats_newtonrap$AIC), 
            (fitstats_henderson$BIC < fitstats_newtonrap$BIC)
        ))
    }
    ### Extract BLUPs
    if (fitstats_comparison >= 0.5) {
        algorithm = "Henderson"
        fitstats = fitstats_henderson
        model = mod_henderson$args
        u = mod_henderson$u[,1]
    } else {
        algorithm = "Newton-Raphson"
        fitstats = fitstats_newtonrap
        model = mod_newtonrap$call
        u = unlist(mod_newtonrap$U)
    }
    loglik = fitstats$logLik
    AIC = fitstats$AIC
    BIC = fitstats$BIC
    return(list(u=u, loglik=loglik, AIC=AIC, BIC=BIC, algorithm=algorithm, model=model))
}
