#' @title rerun_algorithm
#' @description Rerun the algorithm and get the output with the best parameters
#' @param sum_stats dataframe including PRS weights from disease GWAS.
#' @param ped dataframe including phenotype and covariates.
#' @param G The genotype matrix of individual x SNP. Row names are individual ID and column names are SNP ID. 
#' @param bim Allele info for the SNPs in genotype file. The following columns are required: SNP, Ref. Column names are needed.
#' @param initial Indicating the starting values of the predictive effects, 'PRS' or 'zero'.
#' @param num_snp number of total snps with non-zero effects. 
#' It is used to determine the learning rate as 
#' (factor * 1 / num_snp), where the factor is varied across a grid of values: 1, 10, 50, 70. 
#' @param best_criterion The criterion to maximize overall R2 ("best_R2") or the conditional R2 of GxT (conditional on G) ("best_condR2") when selecting the best parameters. Default is "best_R2".
#' @param s2_results The output from the "parameter tuning". 
#' @param covar covariates to be adjusted in the model (other than Tr), such as "age,gender". This is optional. Default is NULL. 
#' @param fixG logical, indicating whether the genetic main (G) effect is fixed. Default is FALSE.
#' @return A list of beta values and test scoresum values.
#' @export
#' @examples
#' \dontrun{
#' data(sum_stats)
#' data(ped)
#' data(G)
#' data(bim)
#' data(parameter_tuning_output)
#' result_list = rerun_algorithm(sum_stats = sum_stats, ped = ped, G = G, bim = bim, 
#' initial = "PRS", num_snp = 2349, s2_results = parameter_tuning_output)
#' }

rerun_algorithm <- function(sum_stats, ped, G, bim, initial, num_snp, s2_results, best_criterion = "best_R2", covar = NULL, fixG = FALSE) {

    if (length(covar) == 0) {
        covar <- 'Tr'
    } else {
        covars <- unlist(strsplit(covar, ","))
        covar <- paste0(c('Tr', covars), collapse = " + ")
    }

    ped$prs_gt <- ped$Tr * ped$prs_g
    G2 <- G  # all samples
    G <- G2[rownames(G2) %in% ped$ID, ]  # full-training samples
    
    # Match SNP and flip +/-
    bim2 <- bim[, c('SNP', 'Ref')]
    colnames(bim2)[2] <- 'Ref_G'
    bim_sum_stats <- inner_join(sum_stats, bim2, by = "SNP")
    bim_sum_stats$Beta2 <- NA
    
    flag1 <- which(bim_sum_stats$Ref == bim_sum_stats$Ref_G)
    if (length(flag1) > 0) bim_sum_stats$Beta2[flag1] <- bim_sum_stats$Beta[flag1]
    
    flag2 <- which(bim_sum_stats$Eff == bim_sum_stats$Ref_G)
    if (length(flag2) > 0) bim_sum_stats$Beta2[flag2] <- -bim_sum_stats$Beta[flag2]
    
    bim_sum_stats <- bim_sum_stats[!is.na(bim_sum_stats$Beta2),]
    
    # Keep the best beta and save
    if (fixG) {
        best <- s2_results[["best_R2"]]
    } else {
        best <- s2_results[[best_criterion]]
    }
    
    best_lambda <- best[1, 'lambda']
    best_factor <- best[1, 'factor']
    best_step <- best[1, 'steps']
    obj <- adjust_y(ped, NULL, covar, fixG)  
    GD_obj <- GD(obj, G, bim_sum_stats, if (fixG) 200 else 30, initial, num_snp, best_lambda, best_factor, fixG)  
    if (best_step == 0) {
        beta <- GD_obj$oc[, c(1:4, 6)]
    } else {
        beta <- GD_obj$oc[, c(1:4, 6 + best_step - 1)]
    }
    colnames(beta)[5] <- 'best_coef'
    
    if (fixG) {
        oc_gt <- beta
        collect_beta <- list(beta_gt = oc_gt)
        #save(collect_beta, file = paste0(s3_output, '_Beta_initial_', initial, '.Rdata'))
    } else {
        oc_g <- beta[1:(nrow(beta) / 2), ]
        oc_gt <- beta[(nrow(beta) / 2 + 1):nrow(beta),]
        collect_beta <- list(beta_g = oc_g, beta_gt = oc_gt)
        #save(collect_beta, file = paste0(s3_output, '_Beta_initial_', initial, '_', best_criterion, '.Rdata'))
    }
    
    # Optional: Apply the beta on the testing to get the scoresum
    if (fixG) {
        snps <- sub("_gt$", "", rownames(oc_gt))
    } else {
        snps <- rownames(oc_g)
    }
    
    geno <- G2[!rownames(G2) %in% ped$ID, snps]  # Keep testing samples not used in full-training set
    
    if (nrow(geno) > 0) {
        # Impute NA as mean
        if (sum(is.na(geno)) > 0) {
            for (j in 1:ncol(geno)) {
                flag <- which(is.na(geno[, j]))
                if (length(flag) > 0) geno[flag, j] <- mean(geno[, j], na.rm = TRUE)
            }
        }
        
        geno <- as.matrix(geno)
        
        if (fixG) {
            weight_gt <- oc_gt$best_coef / oc_gt$sd
            score_gt <- geno %*% weight_gt
            test <- list(score_gt = score_gt)
            #save(test, file = paste0(s3_output, '_testScoresum_initial_', initial, '.Rdata'))
        } else {
            weight_g <- oc_g$best_coef / oc_g$sd
            weight_gt <- oc_gt$best_coef / oc_gt$sd
            score_g <- geno %*% weight_g
            score_gt <- geno %*% weight_gt
            test <- list(score_g = score_g, score_gt = score_gt)
            #save(test, file = paste0(s3_output, '_testScoresum_initial_', initial, '_', best_criterion, '.Rdata'))
        }
    }
    return(list("beta" = collect_beta, "testScoresum" = test))
}


