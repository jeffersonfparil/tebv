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
fn_henderson_vs_newtonraphson_fit = function(mod_henderson, mod_newtonrap, extract_BLUPs=TRUE, verbose=FALSE) {
    ### TEST #################################################################################
    # G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=2, verbose=TRUE)
    # list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=2, dist_effects="norm", n_effects=100, h2=0.75, pheno_reps=3, verbose=TRUE)
    # Y = list_Y_b_E_b_epi$Y
    # df = data.frame(y=as.vector(Y), id=rep(rownames(G), times=3))
    # mod_henderson = tryCatch(sommer::mmec(y ~ 1, random = ~ id, data=df, dateWarning=FALSE, verbose=verbose), error=function(e){NA})
    # mod_newtonrap = tryCatch(sommer::mmer(y ~ 1, random = ~ id, data=df, dateWarning=FALSE, verbose=verbose), error=function(e){NA})
    # extract_BLUPs=TRUE; verbose=TRUE
    ### TEST #################################################################################
    if (is.na(mod_henderson[1]) & is.na(mod_newtonrap[1])) {
        print("Error: unable to fit any model.")
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
    ### Extract BLUPs or BLUEs
    if (fitstats_comparison >= 0.5) {
        ##########################
        ### HENDERSON'S MODELS ###
        ##########################
        algorithm = "Henderson"
        fitstats = fitstats_henderson
        model = mod_henderson$args
        fit = mod_henderson
        ### Sanity check
        if ((sum(grepl("id", model)) == 0) | (sum(grepl("y", model$fixed)) == 0)) {
            print("Please use conventional variable names in the input data frame, i.e. y for the response variables and id for the entry or genotype names.")
            return(NA)
        }
        vec_id = as.character(unique(mod_henderson$data$id))
        vec_environ = as.character(unique(mod_henderson$data$environ))
        vec_environ_x_id = as.character(unique(paste0(mod_henderson$data$environ, ":", mod_henderson$data$id)))
        n_ids = length(vec_id)
        n_environs = length(vec_environ)
        if (extract_BLUPs) {
            #####################
            ### EXTRACT BLUPS ###
            #####################
            u = mod_henderson$u[,1]
            idx = which(names(u) %in% c(vec_id, vec_environ_x_id))
            u = u[idx]
            if (!grepl("rrc\\(environ, id, y, nPC", paste(model$random, collapse="")) & grepl("vsc\\(usc\\(environ\\), isc\\(id\\)", paste(model$random, collapse=""))) {
                if (length(u) == (n_ids*n_environs)) {
                    names(u) = paste0(rep(unique(mod_henderson$data$environ), each=n_ids), ":", names(u))
                } else  if (length(u) == (n_ids + (n_ids*n_environs))) {
                    names(u) = paste0(rep(c("", as.character(unique(mod_henderson$data$environ))), each=n_ids), ":", names(u))
                    names(u) = gsub("^:", "", names(u))
                } else {
                    print("Error: unaccounted for number of effects.")
                    return(NA)
                }
            } else if (grepl("rrc\\(environ, id, y, nPC", paste(model$random, collapse=""))) {
                n_PCs = as.numeric(gsub(" ", "", unlist(strsplit(unlist(strsplit(paste(model$random, collapse=""), "nPC ="))[2], ")"))[1]))
                if (length(u) == (n_ids*n_PCs)) {
                    names(u) = paste0(rep(paste0("PC", 1:n_PCs), each=n_ids), ":", names(u))
                }
                if (length(u) == (n_ids + (n_ids*n_PCs))) {
                    names(u) = paste0(rep(c("", as.character(paste0("PC", 1:n_PCs))), each=n_ids), ":", names(u))
                    names(u) = gsub("^:", "", names(u))
                }
            }
            if (length(u)==0) {
                return(NA)
            }
            if (is.null(u[1])) {
                return(NA)
            }
            V = mod_henderson$sigma ### named vector or matrix for some reason.... (TODO get to the bottom of this!)
        } else {
            #####################
            ### EXTRACT BLUES ###
            #####################
            ### NOTE: BLUEs will never be used for multi-environmental trial analyses
            ### Add the intercept
            intercept = mod_henderson$b[1,1]
            b = mod_henderson$b[,1] + intercept
            ### Define the intercept as the first level of the fixed effects
            b[1] = b[1] - intercept
            control_name = levels(mod_henderson$data$id)[1]
            if (is.null(control_name)) {
                control_name = unique(mod_henderson$data$id)[which(!(unique(mod_henderson$data$id) %in% names(b)))]
            }
            names(b)[1] = control_name
            ### Extract the variance-covariance matrix of BLUEs
            V = mod_henderson$Ci[1:length(b), 1:length(b)]
            rownames(V) = colnames(V) = names(b)
        }
    } else {
        #############################
        ### NEWTON-RAPHSON MODELS ###
        #############################
        algorithm = "Newton-Raphson"
        fitstats = fitstats_newtonrap
        model = mod_newtonrap$call
        fit = mod_newtonrap
        ### Sanity check
        if ((sum(grepl("id", model)) == 0) | (sum(grepl("y", model$fixed)) == 0)) {
            print("Please use conventional variable names in the input data frame, i.e. y for the response variables and id for the entry or genotype names.")
            return(NA)
        }
        n_ids = length(unique(mod_newtonrap$data$id))
        n_environs = length(unique(mod_newtonrap$data$environ))
        if (extract_BLUPs) {
            #####################
            ### EXTRACT BLUPS ###
            #####################
            if (!grepl("rrc\\(environ, id, y, nPC", paste(model$random, collapse=""))) {
                if (grepl("environ:id", paste(model$fixed, collapse="")) | grepl("environ:id", paste(model$random, collapse=""))) {
                    u = unlist(mod_newtonrap$U$`environ:id`$y)
                    names(u) = gsub("^environ", "", names(u))
                    names(u) = gsub(":id", ":", names(u))
                    V = mod_newtonrap$VarU$`environ:id`$y
                } else if (grepl("vsr\\(usr\\(environ\\), id\\)", paste(model$random, collapse=""))) {
                    idx = intersect(grep("id$", names(mod_newtonrap$U)), grep("[^:]*:[^:]*:", names(mod_newtonrap$U), invert=TRUE))
                    u = c()
                    for (i in idx) {
                        u_tmp = mod_newtonrap$U[[i]]$y
                        names(u_tmp) = paste0(gsub("id$", "", names(mod_newtonrap$U)[i]), names(u_tmp))
                        u = c(u, u_tmp)
                    }
                    V = matrix(0, nrow=length(u), ncol=length(u))
                    counter = 1
                    for (i in idx) {
                        V_tmp = mod_newtonrap$VarU[[i]]$y
                        idx_ini = (counter - 1) * nrow(V_tmp) + 1
                        idx_fin = (counter - 0) * nrow(V_tmp)
                        V[idx_ini:idx_fin, idx_ini:idx_fin] = V_tmp
                    }
                } else {
                    u = unlist(mod_newtonrap$U$id$y)
                    if (length(u)==0) {
                        return(NA)
                    }
                    if (is.null(u[1])) {
                        return(NA)
                    }
                    names(u) = gsub("^id", "", names(u))
                    V = mod_newtonrap$VarU$id$y
                }
            } else {
                n_PCs = as.numeric(gsub(" ", "", unlist(strsplit(unlist(strsplit(paste(model$random, collapse=""), "nPC ="))[2], ")"))[1]))
                u = c()
                for (pc in 1:n_PCs) {
                    u = c(u, eval(parse(text=paste0("unlist(mod_newtonrap$U$`PC", pc, ":id`$y)"))))
                }
                if (length(u) == (n_ids*n_PCs)) {
                    names(u) = paste0(rep(paste0("PC", 1:n_PCs), each=n_ids), ":", names(u))
                }
                if (length(u) == (n_ids + (n_ids*n_PCs))) {
                    names(u) = paste0(rep(c("", as.character(paste0("PC", 1:n_PCs))), each=n_ids), ":", names(u))
                    names(u) = gsub("^:", "", names(u))
                }
                V = matrix(0, nrow=length(u), ncol=length(u))
                counter = 1
                for (pc in 1:n_PCs) {
                    V_tmp = eval(parse(text=paste0("unlist(mod_newtonrap$VarU$`PC", pc, ":id`$y)")))
                    idx_ini = (counter - 1) * nrow(V_tmp) + 1
                    idx_fin = (counter - 0) * nrow(V_tmp)
                    V[idx_ini:idx_fin, idx_ini:idx_fin] = V_tmp
                }
            }
            if (is.null(u)) {
                print("No random effects corresponding to the entries found.")
                print("Please consider using `extract_BLUPs=FALSE`.")
                return(NA)
            }
            rownames(V) = colnames(V) = names(u)
        } else {
            #####################
            ### EXTRACT BLUES ###
            #####################
            ### NOTE: BLUEs will never be used for multi-environmental trial analyses
            ### Add the intercept
            intercept = mod_newtonrap$Beta$Estimate[1]
            b = mod_newtonrap$Beta$Estimate + intercept
            names(b) = gsub("^id", "", mod_newtonrap$Beta$Effect)
            ### Define the intercept as the first level of the fixed effects
            b[1] = b[1] - intercept
            control_name = levels(mod_newtonrap$data$id)[1]
            if (is.null(control_name)) {
                control_name = unique(mod_newtonrap$data$id)[which(!(unique(mod_newtonrap$data$id) %in% names(b)))]
            }
            names(b)[1] = control_name
            ### Extract the variance-covariance matrix of BLUEs
            V = mod_newtonrap$VarBeta
            rownames(V) = colnames(V) = names(b)
        }
    }
    loglik = fitstats$logLik
    AIC = fitstats$AIC
    BIC = fitstats$BIC
    if (extract_BLUPs) {
        return(list(u=u, V=V, loglik=loglik, AIC=AIC, BIC=BIC, algorithm=algorithm, model=model, fit=fit))
    } else {
        return(list(b=b, V=V, loglik=loglik, AIC=AIC, BIC=BIC, algorithm=algorithm, model=model, fit=fit))
    }
}

### Extract data frame of GxE breeding values
fn_extract_gxe_breeding_values = function(list_u_V_fitstats) {
    if (is.na(list_u_V_fitstats[1])) {
        print("Error: unable to fit any model.")
        return(NA)
    }
    vec_id = unique(list_u_V_fitstats$fit$data$id)
    vec_environ = unique(list_u_V_fitstats$fit$data$environ)
    vec_environ_x_id = unique(paste0(list_u_V_fitstats$fit$data$environ, ":", list_u_V_fitstats$fit$data$id))
    n = length(vec_id)
    m = length(vec_environ)
    u = list_u_V_fitstats$u
    if (!grepl("rrc\\(environ, id, y, nPC", paste(list_u_V_fitstats$model$random, collapse=""))) {
        ### NOT FACTOR ANALYTIC
        ### Retain only the id and environ:id effects
        idx = which(names(u) %in% c(vec_id, vec_environ_x_id))
        u = u[idx]
        if (length(u) == (n*m)) {
            ### ~ environ:id
            u = u
        } else if (length(u) == (n + (n*m))) {
            ### ~ id + environ:id
            u_id = u[grep(":", names(u), invert=TRUE)]
            u = u[grep(":", names(u), invert=FALSE)]
            for (i in 1:length(u)) {
                # i = 1
                idx = grep(unlist(strsplit(names(u[i]), ":"))[2], names(u_id))
                u[i] = u_gxe[i] + u_id[idx]
            }
        } else {
            print("Error: unaccounted for number of effects.")
            return(NA)
        }
    } else {
        ### FACTOR ANALYTIC
        n_PC = as.numeric(gsub(" ", "", unlist(strsplit(unlist(strsplit(paste(list_u_V_fitstats$fit$args$random, collapse=""), "nPC ="))[2], ")"))[1]))
        if (is.na(n_PC)) {
            n_PC = as.numeric(gsub(" ", "", unlist(strsplit(unlist(strsplit(paste(list_u_V_fitstats$fit$call$random, collapse=""), "nPC ="))[2], ")"))[1]))
        }
        vec_pc_x_id = c()
        for (i in 1:n_PC) {
            vec_pc_x_id = c(vec_pc_x_id, paste0("PC", i, ":", vec_id))
        }
        idx = which(names(u) %in% c(vec_id, vec_pc_x_id))
        u = u[idx]
        vec_pc = unlist(lapply(strsplit(names(u), ":"), FUN=function(x){x[1]}))
        vec_id = unlist(lapply(strsplit(names(u), ":"), FUN=function(x){x[2]}))
        n_PC = length(unique(vec_pc))
        FA_gamma = with(list_u_V_fitstats$fit$data, sommer::rrc(timevar=environ, idvar=id, response=y, nPC=n_PC, returnGamma = TRUE))$Gamma
        df_scores = data.frame(id=unique(vec_id))
        for (i in 1:n_PC) {
            # i = 1
            idx = which(vec_pc == vec_pc[i])
            df_scores = merge(df_scores, data.frame(id=vec_id[idx], pc=u[idx]), by="id")
        }
        mat_scores = as.matrix(df_scores[, 2:(n_PC+1)]); rownames(mat_scores) = df_scores$id
        BLUPs = mat_scores %*% t(FA_gamma)
        u = c()
        u_names = c()
        for (j in 1:ncol(BLUPs)) {
            for (i in 1:nrow(BLUPs)) {
                u = c(u, BLUPs[i, j])
                u_names = c(u_names, paste0(colnames(BLUPs)[j], ":", rownames(BLUPs)[i]))
            }
        }
        names(u) = u_names
    }
    vec_env = unlist(lapply(strsplit(names(u), ":"), FUN=function(x){x[1]}))
    vec_id = unlist(lapply(strsplit(names(u), ":"), FUN=function(x){x[2]}))
    df_BLUPs_GXE = data.frame(id=vec_id, env=vec_env, y=u)
    rownames(df_BLUPs_GXE) = NULL
    return(df_BLUPs_GXE)
}
