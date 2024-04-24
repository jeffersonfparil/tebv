# devtools::load_all()
# library(testthat)

test_that(
    "fn_assess_GXE", {
        print("fn_assess_GXE:")
        ### Simulate data
        set.seed(123)
        n = 100
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=1.0, frac_additional_QTL_per_env=0.15, n_reps=5, verbose=FALSE)
        ### Assess GXE
        list_corr_pval = fn_assess_GXE(list_df_CORR$df, trait="y", id="gen", env="env", verbose=FALSE)
        ### We expect to get similar correlation matrix as the one generated using the simulation library
        x = as.vector(list_corr_pval$mat_corr)
        y = as.vector(list_df_CORR$CORR)
        idx = which(x != 1.0)
        txtplot::txtplot(x[idx], y[idx])
        expect_equal((cor(x[idx], y[idx]) > 0.0), TRUE)
    }
)

test_that(
    "fn_GXE_CRD_BLUPs", {
        print("fn_GXE_CRD_BLUPs:")
        ### Simulate data
        set.seed(456)
        n = 100
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3), env_factor_effects_sd=1.0, frac_additional_QTL_per_env=0.1, n_reps=5, verbose=FALSE)
        # list_df_CORR_list_Y = simquantgen::fn_simulate_gxe(G=G, dist_effects="norm", n_effects=50, purely_additive=TRUE, h2=0.75, env_factor_levels=c(2, 3, 2, 2, 3), env_factor_effects_sd=1.0, frac_additional_QTL_per_env=0.1, n_reps=1, verbose=FALSE)
        # df = list_df_CORR_list_Y$df; trait="y"; id="gen"; env="env"; verbose=TRUE
        # mat_design_matrix_crd_gxe = stats::model.matrix(y ~ gen:env, data=df)
        # dim(mat_design_matrix_crd_gxe)
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_df_CORR_list_Y$list_Y_b_E_b_epi$b))

        # saveRDS(list_df_CORR_list_Y, file="/group/pasture/Jeff/tebv/test_fn_GXE_CRD_BLUPs_list_df_CORR_list_Y.Rds")
        # saveRDS(df_true, file="/group/pasture/Jeff/tebv/test_fn_GXE_CRD_BLUPs_df_true.Rds")

        # list_df_CORR_list_Y = readRDS("/group/pasture/Jeff/tebv/test_fn_GXE_CRD_BLUPs_list_df_CORR_list_Y.Rds")


        ### Fit completely randomised design in a single environment            
        list_CRD_BLUPs = fn_GXE_CRD_BLUPs(df=list_df_CORR_list_Y$df, trait="y", id="gen", env="env", skip_algo_based_on_dimensions=TRUE, verbose=FALSE)
        print(list_CRD_BLUPs)
        ### We expect the true genotype values match with the BLUPs in at least one environment
        merged_BLUPs = merge(list_CRD_BLUPs$df_BLUPs, df_true, by="gen")
        x = merged_BLUPs$y_true
        txtplot::txtdensity(x)
        bool_at_least_one_positive_correlation = FALSE
        for (j in 2:(ncol(merged_BLUPs)-1)) {
            y = merged_BLUPs[,j]
            txtplot::txtdensity(y)
            corr = stats::cor.test(x, y)
            txtplot::txtplot(x, y, xlab="obs", ylab="pred")
            print(paste0("corr = ", round(corr$estimate*100), "% (pval=", round(corr$p.value, 4), ")"))
            if (corr$estimate > 0.0) {
                bool_at_least_one_positive_correlation = TRUE
            }
        }
        expect_equal(bool_at_least_one_positive_correlation, TRUE)
        ### We also expect to capture similar correlations in the BLUPs and the simulated phenotype data
        x = as.vector(list_df_CORR_list_Y$CORR)
        y = as.vector(cor(as.matrix(list_CRD_BLUPs$df_BLUPs[, -1])))
        idx = which(abs(x) < 0.99)
        txtplot::txtplot(x[idx], y[idx])
        expect_equal((cor(x[idx], y[idx]) > 0.0), TRUE)
    }
)