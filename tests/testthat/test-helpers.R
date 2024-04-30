# devtools::load_all()
# library(sommer)
# library(testthat)

test_that(
    "fn_divide_into_lower_and_upper", {
        print("fn_divide_into_lower_and_upper:")
        ### Simulate data
        G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=5, verbose=FALSE)
        # txtplot::txtdensity(list_df_CORR$df$y)
        ### Create a bimodal phenotype distribution by adding 5 to the phenotype values of 20% of the genotypes
        vec_gen_add_10 = sort(sample(unique(list_df_CORR$df$gen), size=20, replace=FALSE))
        idx = which(list_df_CORR$df$gen %in% vec_gen_add_10)
        list_df_CORR$df$y[idx] = list_df_CORR$df$y[idx] + 10.0
        # txtplot::txtdensity(list_df_CORR$df$y)
        ### Divide
        df = fn_divide_into_lower_and_upper(df=list_df_CORR$df, trait="y", threshold=5, verbose=FALSE)
        expect_equal(sum(!is.na(df$y_LOWER)), 2400)
        expect_equal(sum(!is.na(df$y_UPPER)), 600)
    }
)

test_that(
    "fn_find_best_fit_within_algo", {
        print("fn_find_best_fit_within_algo:")
        ### Simulate data
        set.seed(123)
        n_alleles = 2
        pheno_reps = 3
        G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=n_alleles, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, h2=0.75, pheno_reps=pheno_reps, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df = data.frame(y=as.vector(Y), id=rep(rownames(G), times=pheno_reps))
        ### True genotype values
        df_true = data.frame(id=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit using mmec and mmer
        mod_henderson_1 = tryCatch(suppressWarnings(sommer::mmec(y ~ 1, random = ~ id, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        mod_henderson_2 = tryCatch(suppressWarnings(sommer::mmec(y ~ 0, random = ~ id, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        mod_newtonrap_1 = tryCatch(suppressWarnings(sommer::mmer(y ~ 1, random = ~ id, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        mod_newtonrap_2 = tryCatch(suppressWarnings(sommer::mmer(y ~ 0, random = ~ id, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        ### Compare and extract BLUPs
        list_best_mods = fn_find_best_fit_within_algo(
            list_mod_henderson=list(mod_henderson_1, mod_henderson_2), 
            list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2))
        expect_equal(paste(as.character(list_best_mods$mod_henderson$args$random), collapse=""), "~id")
        expect_equal(paste(as.character(list_best_mods$mod_newtonrap$call$random), collapse=""), "~id")
    }
)

test_that(
    "fn_henderson_vs_newtonraphson_fit", {
        print("fn_henderson_vs_newtonraphson_fit:")
        ###########################################
        ### SINGLE ENVIRONMENT
        ###########################################
        ### Simulate data
        set.seed(123)
        n = 100
        n_reps = 3
        n_rows = 10
        n_cols = ceiling(n*n_reps / n_rows)
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=2, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=n_reps, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df = data.frame(y=as.vector(Y), id=rep(rownames(G), times=n_reps), expand.grid(row=1:n_rows, col=1:n_cols))
        ### Simulate row and column effects
        vec_row_effects = rnorm(n=n_rows)
        vec_col_effects = rnorm(n=n_cols)
        for (i in 1:n_rows) {
            for (j in 1:n_cols) {
                idx = which((df$row==i) & (df$col==j))
                df$y[idx] =  df$y[idx] + vec_row_effects[i] + vec_col_effects[j]
            }
        }
        df$row_factor = as.factor(df$row)
        df$col_factor = as.factor(df$col)
        df$block = df$row
        verbose = FALSE
        ### True genotype values
        df_true = data.frame(id=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### List of various models
        vec_mod = suppressWarnings(list(
            mod_henderson_1 = sommer::mmec(y ~ 1, random = ~ id, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_2 = sommer::mmec(y ~ 1 + id, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_3 = sommer::mmec(y ~ 1 + block, random = ~ id, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_4 = sommer::mmec(y ~ 1, random = ~ id + block, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_5 = sommer::mmec(y ~ 1 + id, random = ~ block, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_6 = sommer::mmec(y ~ 1 + id + block, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_7 = sommer::mmec(y ~ 1, random = ~id + row + col + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_8 = sommer::mmec(y ~ 1, random = ~id + row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_9 = sommer::mmec(y ~ 1 + row + col, random = ~id + row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_10 = sommer::mmec(y ~ 1 + row + col + row:col, random = ~id + row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_11 = sommer::mmec(y ~ 1 + row_factor + col_factor + row_factor:col_factor, random = ~id + row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_12 = sommer::mmec(y ~ 1 + id, random = ~row + col + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_13 = sommer::mmec(y ~ 1 + id, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_14 = sommer::mmec(y ~ 1 + id + row + col, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_15 = sommer::mmec(y ~ 1 + id + row + col + row:col, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_henderson_16 = sommer::mmec(y ~ 1 + id + row_factor + col_factor + row_factor:col_factor, random = ~row_factor + col_factor + row_factor:col_factor, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_1 = sommer::mmer(y ~ 1, random = ~ id, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_2 = sommer::mmer(y ~ 1 + id, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_3 = sommer::mmer(y ~ 1 + block, random = ~ id, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_4 = sommer::mmer(y ~ 1, random = ~ id + block, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_5 = sommer::mmer(y ~ 1 + id, random = ~ block, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_6 = sommer::mmer(y ~ 1 + id + block, rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_7 = sommer::mmer(y ~ 1, random = ~id + row + col + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_8 = sommer::mmer(y ~ 1, random = ~id + row_factor + col_factor + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_9 = sommer::mmer(y ~ 1 + row_factor + col_factor, random = ~id + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_10 = sommer::mmer(y ~ 1 + row + col + row:col, random = ~id + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_11 = sommer::mmer(y ~ 1 + row + col + row:col, random = ~id + row_factor + col_factor + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_12 = sommer::mmer(y ~ 1 + id, random = ~row + col + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_13 = sommer::mmer(y ~ 1 + id, random = ~row_factor + col_factor + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_14 = sommer::mmer(y ~ 1 + id + row_factor + col_factor, random = ~sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_15 = sommer::mmer(y ~ 1 + id + row + col + row:col, random = ~sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose),
            mod_newtonrap_16 = sommer::mmer(y ~ 1 + id + row + col + row:col, random = ~row_factor + col_factor + sommer::spl2Da(x.coord=df$col, y.coord=df$row, nsegments=c(n_cols, n_rows), degree=c(3,3)), rcov= ~ units, data=df, dateWarning=FALSE, verbose=verbose)
        ))
        ### Assess fit
        for (i in 1:length(vec_mod)) {
            mod = vec_mod[[i]]
            algo = names(vec_mod)[i]
            print("###################################")
            print(paste0(i, " | ", algo))
            if (grepl("henderson", algo)) {
                list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=mod, mod_newtonrap=NA, extract_BLUPs=TRUE, verbose=verbose))
                if ((is.na(list_ub_V_fitstats[1])) | (tryCatch(sum(grepl("entry", names(list_ub_V_fitstats$u)))==0, error=function(e){TRUE}))) {
                    list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=mod, mod_newtonrap=NA, extract_BLUPs=FALSE, verbose=verbose))
                    df_effects = data.frame(id=names(list_ub_V_fitstats$b), y=list_ub_V_fitstats$b); rownames(df_effects) = NULL
                } else {
                    df_effects = data.frame(id=names(list_ub_V_fitstats$u), y=list_ub_V_fitstats$u); rownames(df_effects) = NULL
                }
            } else {
                list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod, extract_BLUPs=TRUE, verbose=verbose))
                if ((is.na(list_ub_V_fitstats[1])) | (tryCatch(sum(grepl("entry", names(list_ub_V_fitstats$u)))==0, error=function(e){TRUE}))) {
                    list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod, extract_BLUPs=FALSE, verbose=verbose))
                    df_effects = data.frame(id=names(list_ub_V_fitstats$b), y=list_ub_V_fitstats$b); rownames(df_effects) = NULL
                } else {
                    df_effects = data.frame(id=names(list_ub_V_fitstats$u), y=list_ub_V_fitstats$u); rownames(df_effects) = NULL
                }
            }
            merged_effects = merge(df_effects, df_true, by="id")
            corr = cor(merged_effects$y, merged_effects$y_true)
            txtplot::txtplot(merged_effects$y_true, merged_effects$y, xlab="obs", ylab="pred")
            print(paste0("corr=", round(100*corr), "%"))
            expect_equal(corr > 0, TRUE)
        }



        ###########################################
        ### MULTIPLE ENVIRONMENTS
        ###########################################
        ### Simulate data
        set.seed(123)
        n = 100
        n_reps = 5
        n_rows = 10
        n_cols = ceiling(n*n_reps / n_rows)
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=n_reps, verbose=FALSE)
        df = list_df_CORR_list_Y$df
        df_row_col = expand.grid(row=1:n_rows, col=1:n_cols)
        df$row = rep(df_row_col$row, times=length(unique(df$env)))
        df$col = rep(df_row_col$col, times=length(unique(df$env)))
        ### Simulate row and column effects
        for (env in unique(df$env)) {
            vec_row_effects = rnorm(n=n_rows)
            vec_col_effects = rnorm(n=n_cols)
            for (i in 1:n_rows) {
                for (j in 1:n_cols) {
                    idx = which((df$env == env) & (df$row==i) & (df$col==j))
                    df$y[idx] =  df$y[idx] + vec_row_effects[i] + vec_col_effects[j]
                }
            }
        }
        df = data.frame(
            y = df$y, 
            id = df$gen, 
            environ = df$env,
            row = df$row,
            col = df$col
        )
        df$row_factor = as.factor(df$row)
        df$col_factor = as.factor(df$col)
        n_rows = nlevels(df$row_factor)
        n_cols = nlevels(df$col_factor)
        df$block = df$row
        tolParInv = 0.01
        verbose = FALSE

        mod_newtonrap_1 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ environ:id, rcov= ~ units, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_2 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + environ:id, rcov= ~ units, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_3 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ vsr(usr(environ), id), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_4 = tryCatch(sommer::mmer(y ~ 1 + environ, random= ~ id + vsr(usr(environ), id), rcov= ~ vsr(dsr(environ), units), tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_5 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=2)), id), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_newtonrap_6 = tryCatch(sommer::mmer(y ~ 1 + environ, random = ~ vsr(usr(rrc(environ, id, y, nPC=1)), id), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        ### Henderson's models
        mod_henderson_1 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ environ:id, rcov= ~ units, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_2 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + environ:id, rcov= ~ units, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_3 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ vsc(usc(environ), isc(id)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_4 = tryCatch(sommer::mmec(y ~ 1 + environ, random= ~ id + vsc(usc(environ), isc(id)), rcov= ~ vsc(dsc(environ), isc(units)), tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_5 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=2)), isc(id)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})
        mod_henderson_6 = tryCatch(sommer::mmec(y ~ 1 + environ, random = ~ vsc(usc(rrc(environ, id, y, nPC=1)), isc(id)), rcov= ~ units, nIters=100, tolParInv=tolParInv, data=df, dateWarning=FALSE, verbose=verbose),
            error=function(e){NA})


        ### Assess fit
        for (i in 1:length(vec_mod)) {
            mod = vec_mod[[i]]
            algo = names(vec_mod)[i]
            print("###################################")
            print(paste0(i, " | ", algo))
            if (grepl("henderson", algo)) {
                list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=mod, mod_newtonrap=NA, extract_BLUPs=TRUE, verbose=verbose))
                if ((is.na(list_ub_V_fitstats[1])) | (tryCatch(sum(grepl("entry", names(list_ub_V_fitstats$u)))==0, error=function(e){TRUE}))) {
                    list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=mod, mod_newtonrap=NA, extract_BLUPs=FALSE, verbose=verbose))
                    df_effects = data.frame(id=names(list_ub_V_fitstats$b), y=list_ub_V_fitstats$b); rownames(df_effects) = NULL
                } else {
                    df_effects = data.frame(id=names(list_ub_V_fitstats$u), y=list_ub_V_fitstats$u); rownames(df_effects) = NULL
                }
            } else {
                list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod, extract_BLUPs=TRUE, verbose=verbose))
                if ((is.na(list_ub_V_fitstats[1])) | (tryCatch(sum(grepl("entry", names(list_ub_V_fitstats$u)))==0, error=function(e){TRUE}))) {
                    list_ub_V_fitstats = suppressWarnings(fn_henderson_vs_newtonraphson_fit(mod_henderson=NA, mod_newtonrap=mod, extract_BLUPs=FALSE, verbose=verbose))
                    df_effects = data.frame(id=names(list_ub_V_fitstats$b), y=list_ub_V_fitstats$b); rownames(df_effects) = NULL
                } else {
                    df_effects = data.frame(id=names(list_ub_V_fitstats$u), y=list_ub_V_fitstats$u); rownames(df_effects) = NULL
                }
            }
            merged_effects = merge(df_effects, df_true, by="id")
            corr = cor(merged_effects$y, merged_effects$y_true)
            txtplot::txtplot(merged_effects$y_true, merged_effects$y, xlab="obs", ylab="pred")
            print(paste0("corr=", round(100*corr), "%"))
            expect_equal(corr > 0, TRUE)
        }
    }
)