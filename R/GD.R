
# function to apply Gradient descent and keep all beta
GD <- function(obj, G, bim_sum_stats, step = 15, initial, nsnp, lambda = NULL, factor1 = NULL, fixG = FALSE) {
    y_res = obj$y_res
    geno = G[y_res$ID, bim_sum_stats$SNP]
    geno_info = bim_sum_stats[, c('SNP', 'Ref_G', 'Beta2')]
    
    # Impute NA as mean
    if (sum(is.na(geno)) > 0) {
        for (j in 1:ncol(geno)) {
            flag = which(is.na(geno[, j]))
            if (length(flag) > 0) { 
                geno[flag, j] = mean(geno[, j], na.rm = TRUE) 
            }
        }
    }
    
    geno_info$mean = colMeans(geno, na.rm = TRUE)
    geno_info$sd = apply(geno, 2, sd)
    list1 = which(geno_info$sd == 0)
    if (length(list1) > 0) {
        geno_info = geno_info[-list1, ]
        geno = geno[, -list1]
    }
    
    # Add the geno x T (no NA in T)
    table(rownames(geno) == y_res$ID)
    GT = geno * y_res$Tr 
    geno_info$mean_gt = apply(GT, 2, mean)
    geno_info$sd_gt = apply(GT, 2, sd)
   
    table(colnames(GT) == colnames(G))
    colnames(GT) = paste0(colnames(GT), '_gt')
   
    # Standardized 
    for (j in 1:ncol(geno)) {
        geno[, j] = (geno[, j] - geno_info$mean[j]) / geno_info$sd[j]
        GT[, j] = (GT[, j] - geno_info$mean_gt[j]) / geno_info$sd_gt[j]
    }
    
    if (fixG) {
        X = as.matrix(GT)
    } else {
        X = as.matrix(cbind(geno, GT))
    }
    
    # Check
    stopifnot(all(rownames(X) == y_res$ID))
    stopifnot(all(c(sub("_gt$", "", colnames(X))) == geno_info$SNP))

    
    GG = t(X) %*% X / (nrow(X) - 1)
    gy = t(X) %*% (y_res$res) / (nrow(y_res) - 1)
    if (nrow(geno_info) > 0) {
        if (fixG) {
            if (initial == 'PRS') {
                betatemp_gt = geno_info$Beta2 * geno_info$sd_gt * obj$coeff
            }
            if (initial == 'zero') {
                betatemp_gt = rep(0, nrow(geno_info))
            }
            betatemp = betatemp_gt
        } else {
            betatemp_g = geno_info$Beta2 * geno_info$sd * obj$coeff[1]
            if (initial == 'PRS') {
                betatemp_gt = geno_info$Beta2 * geno_info$sd_gt * obj$coeff[2]
            }
            if (initial == 'zero') {
                betatemp_gt = rep(0, length(betatemp_g))
            }
            betatemp = c(betatemp_g, betatemp_gt)
        }
        
        u0 = gy - GG %*% betatemp
        beta.all = cbind(u0, betatemp)
        
        if (is.null(lambda)) {
            lambda_vec = c(0, 0.5, 0.99)
        } else {
            lambda_vec = c(lambda)
        }
        if (is.null(factor1)) {
            factor1_vec = c(1, 10, 50, 70)
        } else {
            factor1_vec = c(factor1)
        }

        tunning = data.frame(lambda = NA, factor = NA, steps = NA)
        for (lambda in lambda_vec) {
            for (factor1 in factor1_vec) {
                k = 1
                betatemp = beta.all[, 2]
                u0 = beta.all[, 1]
                while (k <= step) {
                    learningrate = 1 / nsnp * factor1
                    if (learningrate > 1) {
                        learningrate = 1
                    }
                    #print(learningrate)
                    beta_old = betatemp
                    betatemp = learningrate * u0 + (1 - lambda) * beta_old
                    u0 = u0 - GG %*% (betatemp - beta_old)
                    
                    beta.all = cbind(beta.all, betatemp)
                    k = k + 1
                    
                    tunning = as.data.frame(rbind(tunning, c(lambda, factor1, k)))
                } 
            }
        }
        
        if (fixG) {
            oc = data.frame(
                SNP = rownames(beta.all), 
                Ref_G = geno_info$Ref_G,
                sd = geno_info$sd_gt, 
                mean = geno_info$mean_gt
            )
        } else {
            oc = data.frame(
                SNP = rownames(beta.all), 
                Ref_G = c(geno_info$Ref_G, geno_info$Ref_G),
                sd = c(geno_info$sd, geno_info$sd_gt), 
                mean = c(geno_info$mean, geno_info$mean_gt)
            )
        }
        
        oc = as.data.frame(cbind(oc, beta.all))
        tunning = tunning[-1, ]
        
        return(list(oc = oc, tunning = tunning))
    }
}

