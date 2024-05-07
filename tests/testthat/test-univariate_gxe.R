# devtools::load_all()
# library(sommer)
# library(testthat)

test_that("fn_assess_GXE",
    {
        print("fn_assess_GXE:")
        ### Simulate data
        set.seed(123)
        n = 30
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=1.0, frac_additional_QTL_per_env=0.15, n_reps=5, verbose=FALSE)
        ### Assess GXE
        list_corr_pval = fn_assess_GXE(list_df_CORR$df, trait="y", id="gen", env="env", verbose=FALSE)
        ### We expect to get similar correlation matrix as the one generated using the simulation library
        x = as.vector(list_corr_pval$mat_corr)
        y = as.vector(list_df_CORR$CORR)
        idx = which(x != 1.0)
        corr = cor(x[idx], y[idx])
        print(paste0("corr=", round(100*corr), "%"))
        txtplot::txtplot(x[idx], y[idx])
        expect_equal((corr > 0.0), TRUE)
    }
)

test_that("fn_GXE_CRD_BLUPs",
    {
        print("fn_GXE_CRD_BLUPs:")
        ### Simulate data
        set.seed(123)
        n = 30
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=1.0, frac_additional_QTL_per_env=0.1, n_reps=5, verbose=FALSE)
        ### True genotype values
        df_true = data.frame(id=rownames(G), y_true=(G %*% list_df_CORR_list_Y$list_Y_b_E_b_epi$b))
        ### Fit completely randomised design in a single environment            
        list_CRD_BLUPs = fn_GXE_CRD_BLUPs(df=list_df_CORR_list_Y$df, trait="y", id="gen", env="env", tolParInv=0.01, skip_algo_based_on_dimensions=TRUE, verbose=FALSE)
        ### We expect the true genotype values match with the BLUPs in at least one environment
        merged_BLUPs = merge(list_CRD_BLUPs$df_effects, df_true, by="id")
        x = merged_BLUPs$y_true
        y = merged_BLUPs$y
        txtplot::txtdensity(x)
        txtplot::txtdensity(y)
        corr = stats::cor.test(x, y)
        print(paste0("corr = ", round(corr$estimate*100), "% (pval=", round(corr$p.value, 4), ")"))
        txtplot::txtplot(x, y, xlab="obs", ylab="pred")
        expect_equal((corr$estimate < 0.0)[["cor"]], FALSE)
    }
)

test_that("fn_GXE_RBD_BLUPs",
    {
        print("fn_GXE_RBD_BLUPs:")
        ### Simulate data
        set.seed(123)
        n = 30
        n_blocks = 5
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=1.0, frac_additional_QTL_per_env=0.1, n_reps=n_blocks, verbose=FALSE)
        ### True genotype values
        df_true = data.frame(id=rownames(G), y_true=(G %*% list_df_CORR_list_Y$list_Y_b_E_b_epi$b))
        ### Simulate block effects
        for (env in unique(list_df_CORR_list_Y$df$env)) {
            for (j in 1:n_blocks) {
                idx = which((list_df_CORR_list_Y$df$env == env) & (list_df_CORR_list_Y$df$rep == as.character(j)))
                list_df_CORR_list_Y$df$y[idx] = list_df_CORR_list_Y$df$y[idx] + stats::rnorm(1)
            }
        }
        ### Fit completely randomised design in a single environment            
        list_RBD_BLUPs = fn_GXE_RBD_BLUPs(df=list_df_CORR_list_Y$df, trait="y", id="gen", env="env", block="rep", tolParInv=0.01, skip_algo_based_on_dimensions=TRUE, verbose=FALSE)
        ### We expect the true genotype values match with the BLUPs in at least one environment
        merged_BLUPs = merge(list_RBD_BLUPs$df_effects, df_true, by="id")
        x = merged_BLUPs$y_true
        y = merged_BLUPs$y
        txtplot::txtdensity(x)
        txtplot::txtdensity(y)
        corr = stats::cor.test(x, y)
        print(paste0("corr = ", round(corr$estimate*100), "% (pval=", round(corr$p.value, 4), ")"))
        txtplot::txtplot(x, y, xlab="obs", ylab="pred")
        expect_equal((corr$estimate < 0.0)[["cor"]], FALSE)
    }
)

test_that("fn_GXE_SPAT_BLUPs",
    {
        print("fn_GXE_SPAT_BLUPs:")
        ### Simulate data
        set.seed(123)
        n = 30
        n_reps = 5
        n_rows = 10
        n_cols = n*n_reps / n_rows
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=n_reps, verbose=FALSE)
        df = list_df_CORR_list_Y$df
        df_row_col = expand.grid(row=1:n_rows, col=1:n_cols)
        df$row = rep(df_row_col$row, times=length(unique(df$env)))
        df$col = rep(df_row_col$col, times=length(unique(df$env)))
        ### True genotype values
        df_true = data.frame(id=rownames(G), y_true=(G %*% list_df_CORR_list_Y$list_Y_b_E_b_epi$b))
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
        ### Fit completely randomised design in a single environment            
        list_SPAT_BLUPs = fn_GXE_SPAT_BLUPs(df=df, trait="y", id="gen", env="env", row="row", col="col", tolParInv=0.01, skip_algo_based_on_dimensions=TRUE, verbose=FALSE)
        ### We expect the true genotype values match with the BLUPs in at least one environment
        merged_BLUPs = merge(list_SPAT_BLUPs$df_effects, df_true, by="id")
        x = merged_BLUPs$y_true
        y = merged_BLUPs$y
        txtplot::txtdensity(x)
        txtplot::txtdensity(y)
        corr = stats::cor.test(x, y)
        print(paste0("corr = ", round(corr$estimate*100), "% (pval=", round(corr$p.value, 4), ")"))
        txtplot::txtplot(x, y, xlab="obs", ylab="pred")
        expect_equal((corr$estimate < 0.0)[["cor"]], FALSE)
    }
)