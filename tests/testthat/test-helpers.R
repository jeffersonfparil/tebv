# devtools::load_all()
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
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=pheno_reps))
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit using mmec and mmer
        mod_henderson_1 = tryCatch(suppressWarnings(sommer::mmec(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        mod_henderson_2 = tryCatch(suppressWarnings(sommer::mmec(y ~ 0, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        mod_newtonrap_1 = tryCatch(suppressWarnings(sommer::mmer(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        mod_newtonrap_2 = tryCatch(suppressWarnings(sommer::mmer(y ~ 0, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE)), error=function(e){NA})
        ### Compare and extract BLUPs
        list_best_mods = fn_find_best_fit_within_algo(
            list_mod_henderson=list(mod_henderson_1, mod_henderson_2), 
            list_mod_newtonrap=list(mod_newtonrap_1, mod_newtonrap_2))
        expect_equal(paste(as.character(list_best_mods$mod_henderson$args$random), collapse=""), "~gen")
        expect_equal(paste(as.character(list_best_mods$mod_newtonrap$call$random), collapse=""), "~gen")
    }
)

test_that(
    "fn_henderson_vs_newtonraphson_fit", {
        print("fn_henderson_vs_newtonraphson_fit:")
        ### Simulate data
        set.seed(123)
        n_alleles = 2
        pheno_reps = 3
        G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=n_alleles, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, h2=0.75, pheno_reps=pheno_reps, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=pheno_reps))
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit using mmec and mmer
        mod_henderson = tryCatch(sommer::mmec(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE), error=function(e){NA})
        mod_newtonrap = tryCatch(sommer::mmer(y ~ 1, random = ~ gen, data=df, dateWarning=FALSE, verbose=FALSE), error=function(e){NA})
        ### Compare and extract BLUPs
        list_u_fitstats = fn_henderson_vs_newtonraphson_fit(mod_henderson=mod_henderson, mod_newtonrap=mod_newtonrap, verbose=FALSE)
        ### We expect the true genotype values match closely with the BLUPs
        names(list_u_fitstats$u) = gsub("gen.y.gen", "", names(list_u_fitstats$u))
        df_BLUPs = data.frame(gen=names(list_u_fitstats$u), y=list_u_fitstats$u); rownames(df_BLUPs) = NULL
        merged_BLUPs = merge(df_BLUPs, df_true, by="gen")
        txtplot::txtplot(merged_BLUPs$y_true, merged_BLUPs$y, xlab="obs", ylab="pred")
        expect_equal(round(cor(merged_BLUPs$y, merged_BLUPs$y_true)), 1)
    }
)