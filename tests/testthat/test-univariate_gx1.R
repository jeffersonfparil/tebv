# devtools::load_all()
# library(sommer)
# library(testthat)

test_that(
    "fn_GX1_CRD_BLUPs", {
        print("fn_GX1_CRD_BLUPs:")
        ### Simulate data
        set.seed(123)
        n_alleles = 2
        pheno_reps = 3
        G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=n_alleles, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=pheno_reps, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=pheno_reps))
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit completely randomised design in a single environment            
        list_CRD_BLUPs = suppressWarnings(fn_GX1_CRD_BLUPs(df=df, trait="y", id="gen", verbose=FALSE))
        ### We expect the true genotype values match closely with the BLUPs
        df_BLUPs = list_CRD_BLUPs$df_effects
        merged_BLUPs = merge(df_BLUPs, df_true, by="gen")
        corr = cor(merged_BLUPs$y, merged_BLUPs$y_true)
        txtplot::txtplot(merged_BLUPs$y_true, merged_BLUPs$y, xlab="obs", ylab="pred")
        print(paste0("corr=", round(100*corr), "%"))
        expect_equal(round(corr), 1)
    }
)

test_that(
    "fn_GX1_CRD_BLUEs", {
        print("fn_GX1_CRD_BLUEs:")
        ### Simulate data
        set.seed(123)
        n_alleles = 2
        pheno_reps = 3
        G = simquantgen::fn_simulate_genotypes(n=100, l=1000, ploidy=42, n_alleles=n_alleles, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=pheno_reps, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=pheno_reps))
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit completely randomised design in a single environment            
        list_CRD_BLUEs = suppressWarnings(fn_GX1_CRD_BLUEs(df=df, control_id=rownames(G)[17], trait="y", id="gen", verbose=FALSE))
        ### We expect the true genotype values match closely with the BLUEs
        df_BLUEs = list_CRD_BLUEs$df_effects
        merged_BLUEs = merge(df_BLUEs, df_true, by="gen")
        corr = cor(merged_BLUEs$y, merged_BLUEs$y_true)
        txtplot::txtplot(merged_BLUEs$y_true, merged_BLUEs$y, xlab="obs", ylab="pred")
        print(paste0("corr=", round(100*corr), "%"))
        expect_equal(round(corr), 1)
    }
)

test_that(
    "fn_GX1_RBD_BLUPs", {
        print("fn_GX1_RBD_BLUPs:")
        ### Simulate data
        set.seed(123)
        n = 100
        n_blocks = 10
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=2, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=n_blocks, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Simulate block effects
        for (j in 1:n_blocks) {
            Y[, j] = Y[, j] + rnorm(1)
        }
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=n_blocks), rep=rep(1:n_blocks, each=nrow(Y)))
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit randomised complete block design in a single environment            
        list_RBD_BLUPs = suppressWarnings(fn_GX1_RBD_BLUPs(df=df, trait="y", id="gen", block="rep", verbose=FALSE))
        ### We expect the true genotype values match closely with the BLUPs
        df_BLUPs = list_RBD_BLUPs$df_effects
        merged_BLUPs = merge(df_BLUPs, df_true, by="gen")
        corr = cor(merged_BLUPs$y, merged_BLUPs$y_true)
        txtplot::txtplot(merged_BLUPs$y_true, merged_BLUPs$y, xlab="obs", ylab="pred")
        print(paste0("corr=", round(100*corr), "%"))
        expect_equal(round(cor(merged_BLUPs$y, merged_BLUPs$y_true)), 1)
    }
)

test_that(
    "fn_GX1_RBD_BLUEs", {
        print("fn_GX1_RBD_BLUEs:")
        ### Simulate data
        set.seed(123)
        n = 100
        n_blocks = 10
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=2, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=n_blocks, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Simulate block effects
        for (j in 1:n_blocks) {
            Y[, j] = Y[, j] + rnorm(1)
        }
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=n_blocks), rep=rep(1:n_blocks, each=nrow(Y)))
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit randomised complete block design in a single environment            
        list_RBD_BLUEs = suppressWarnings(fn_GX1_RBD_BLUEs(df=df, control_id=rownames(G)[17], trait="y", id="gen", block="rep", verbose=FALSE))
        ### We expect the true genotype values match closely with the BLUEs
        df_BLUEs = list_RBD_BLUEs$df_effects
        merged_BLUEs = merge(df_BLUEs, df_true, by="gen")
        corr = cor(merged_BLUEs$y, merged_BLUEs$y_true)
        txtplot::txtplot(merged_BLUEs$y_true, merged_BLUEs$y, xlab="obs", ylab="pred")
        print(paste0("corr=", round(100*corr), "%"))
        expect_equal(round(cor(merged_BLUEs$y, merged_BLUEs$y_true)), 1)
    }
)

test_that(
    "fn_GX1_SPAT_BLUPs", {
        print("fn_GX1_SPAT_BLUPs:")
        ### Simulate data
        set.seed(123)
        n = 100
        n_reps = 3
        n_rows = 10
        n_cols = ceiling(n*n_reps / n_rows)
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=2, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=n_reps, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=n_reps), expand.grid(row=1:n_rows, col=1:n_cols))
        ### Simulate row and column effects
        vec_row_effects = rnorm(n=n_rows)
        vec_col_effects = rnorm(n=n_cols)
        for (i in 1:n_rows) {
            for (j in 1:n_cols) {
                idx = which((df$row==i) & (df$col==j))
                df$y[idx] =  df$y[idx] + vec_row_effects[i] + vec_col_effects[j]
            }
        }
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit randomised complete block design in a single environment            
        list_SPAT_BLUPs = suppressWarnings(fn_GX1_SPAT_BLUPs(df=df, trait="y", id="gen", row="row", col="col", verbose=FALSE))
        ### Print the selected/best spatial model
        print(list_SPAT_BLUPs$model)
        ### We expect the true genotype values match closely with the BLUPs
        df_BLUPs = list_SPAT_BLUPs$df_effects
        merged_BLUPs = merge(df_BLUPs, df_true, by="gen")
        corr = cor(merged_BLUPs$y, merged_BLUPs$y_true)
        txtplot::txtplot(merged_BLUPs$y_true, merged_BLUPs$y, xlab="obs", ylab="pred")
        print(paste0("corr=", round(100*corr), "%"))
        expect_equal(round(cor(merged_BLUPs$y, merged_BLUPs$y_true)), 1)
    }
)

test_that(
    "fn_GX1_SPAT_BLUEs", {
        print("fn_GX1_SPAT_BLUEs:")
        ### Simulate data
        set.seed(123)
        n = 100
        n_reps = 3
        n_rows = 10
        n_cols = ceiling(n*n_reps / n_rows)
        G = simquantgen::fn_simulate_genotypes(n=n, l=1000, ploidy=42, n_alleles=2, verbose=FALSE)
        list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=2, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=n_reps, verbose=FALSE)
        Y = list_Y_b_E_b_epi$Y
        df = data.frame(y=as.vector(Y), gen=rep(rownames(G), times=n_reps), expand.grid(row=1:n_rows, col=1:n_cols))
        ### Simulate row and column effects
        vec_row_effects = rnorm(n=n_rows)
        vec_col_effects = rnorm(n=n_cols)
        for (i in 1:n_rows) {
            for (j in 1:n_cols) {
                idx = which((df$row==i) & (df$col==j))
                df$y[idx] =  df$y[idx] + vec_row_effects[i] + vec_col_effects[j]
            }
        }
        ### True genotype values
        df_true = data.frame(gen=rownames(G), y_true=(G %*% list_Y_b_E_b_epi$b))
        ### Fit randomised complete block design in a single environment            
        list_SPAT_BLUEs = suppressWarnings(fn_GX1_SPAT_BLUEs(df=df, control_id=rownames(G)[17], trait="y", id="gen", row="row", col="col", verbose=FALSE))
        ### We expect the true genotype values match closely with the BLUEs
        df_BLUEs = list_SPAT_BLUEs$df_effects
        merged_BLUEs = merge(df_BLUEs, df_true, by="gen")
        corr = cor(merged_BLUEs$y, merged_BLUEs$y_true)
        txtplot::txtplot(merged_BLUEs$y_true, merged_BLUEs$y, xlab="obs", ylab="pred")
        print(paste0("corr=", round(100*corr), "%"))
        expect_equal(round(cor(merged_BLUEs$y, merged_BLUEs$y_true)), 1)
    }
)