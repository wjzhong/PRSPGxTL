#' @title parameter_tuning
#' @description parameter tuning
#' @param X dataframe including phenotype and covariates.
#' @param s1_results a list of output from the "inner_layer_CV" function.
#' @param s1_files The file name for output from the "inner_layer_CV" function, if loading the result from a ".rds" file. 
#' @param s1_byCHR_files a vector or list of file names for output from the "inner_layer_CV" function for each chromosome, if loading the results from ".rds" files. 
#' @param covar covariates to be adjusted in the model (other than Tr), such as "age,gender". This is optional. Default is NULL. 
#' @param fixG logical, indicating whether the genetic main (G) effect is fixed. Default is FALSE.
#' @param verbose logical, indicating whether to print more details. Default is FALSE.
#' @examples
#' \dontrun{
#' data(ped)
#' data(inner_layer_CV_output)
#' result_list = parameter_tuning(X = ped, s1_results = inner_layer_CV_output)
#' }

parameter_tuning <- function(X, s1_results = NULL, s1_files = NULL, s1_byCHR_files = NULL, fixG = FALSE, covar = NULL, verbose = FALSE) {

    params <- list(s1_results, s1_files, s1_byCHR_files)
    num_provided <- sum(!sapply(params, is.null))
    if (num_provided != 1) {
        stop("Exactly one of 's1_results', 's1_files', or 's1_byCHR_files' must be provided.")
    }

    if (length(covar) == 0) {
        covar <- 'Tr'
    } else {
        covars <- unlist(strsplit(covar, ","))
        covar <- paste0(c('Tr', covars), collapse = " + ")
    }

  if (!is.null(s1_files)) {
    s1_results = readRDS(s1_files)
  }

  collect <- function(j) {
    
    if (!is.null(s1_byCHR_files)) {
      result <- readRDS(s1_byCHR_files[[1]])[[j]]
    } else {
      result <- s1_results[[j]]
    }
    
    if (fixG) {
        prs_gt <- result$score_gt
        prs_gt_train <- result$score_gt_train
    } else {
        prs_g <- result$score_g
        prs_gt <- result$score_gt
        prs_g_train <- result$score_g_train
        prs_gt_train <- result$score_gt_train
    }

    tunning <- result$tunning
    tunning <- rbind(rep(0, ncol(tunning)), tunning)
    #add the prs across all chromosomes
    if (!is.null(s1_byCHR_files)) {
      for (chr in 2:22) {
        result <- readRDS(s1_byCHR_files[[chr]])[[j]]
        if (fixG) {
            prs_gt <- prs_gt + result$score_gt
            prs_gt_train <- prs_gt_train + result$score_gt_train
        } else {
            prs_g <- prs_g + result$score_g
            prs_gt <- prs_gt + result$score_gt
            prs_g_train <- prs_g_train + result$score_g_train
            prs_gt_train <- prs_gt_train + result$score_gt_train
        }
      }
    }

    df <- X[X$ID %in% rownames(prs_gt), ]
    if (fixG) {
      prs_gt <- prs_gt[df$ID, ]
    } else {
      prs_g <- prs_g[df$ID, ]
      prs_gt <- prs_gt[df$ID, ]
    }
    prs_gtt <- prs_gt * df$Tr

    df_train <- X[X$ID %in% rownames(prs_gt_train), ]
    if (fixG) {
      prs_gt_train <- prs_gt_train[df_train$ID, ]
    } else {
      prs_g_train <- prs_g_train[df_train$ID, ]
      prs_gt_train <- prs_gt_train[df_train$ID, ]
    }
    prs_gtt_train <- prs_gt_train * df_train$Tr

    model0 <- lm(paste0('Y ~ ', covar), data = df)
    df$res <- model0$residuals
    model0_train <- lm(paste0('Y ~ ', covar), data = df_train)
    df_train$res <- model0_train$residuals

    tunning$R2 <- tunning$R2_g <- tunning$R2_gt <- tunning$R2_train <- tunning$R2_g_train <- tunning$R2_gt_train <- tunning$condR2_gt <- tunning$condR2_gt_train <- NA
    for (i in 1:nrow(tunning)) {
      if (fixG) {
        df$prs_gtt <- prs_gtt[, i]
        model1 <- lm(res ~ prs_g + prs_gtt, data = df)
        tunning[i, 'R2'] <- summary(model1)$r.squared
        model2 <- lm(res ~ prs_g, data = df)
        tunning[i, 'R2_g'] <- summary(model2)$r.squared
        model3 <- lm(res ~ prs_gtt, data = df)
        tunning[i, 'R2_gt'] <- summary(model3)$r.squared
        model_g <- lm(paste0('Y ~ prs_g + ', covar), data = df)
        df$res_g <- model_g$residuals
        model4 <- lm(res_g ~ prs_gtt, data = df)
        tunning[i, 'condR2_gt'] <- summary(model4)$r.squared

        df_train$prs_gtt <- prs_gtt_train[, i]
        model1 <- lm(res ~ prs_g + prs_gtt, data = df_train)
        tunning[i, 'R2_train'] <- summary(model1)$r.squared
        model2 <- lm(res ~ prs_g, data = df_train)
        tunning[i, 'R2_g_train'] <- summary(model2)$r.squared
        model3 <- lm(res ~ prs_gtt, data = df_train)
        tunning[i, 'R2_gt_train'] <- summary(model3)$r.squared
        model_g <- lm(paste0('Y ~ prs_g + ', covar), data = df_train)
        df_train$res_g <- model_g$residuals
        model4 <- lm(res_g ~ prs_gtt, data = df_train)
        tunning[i, 'condR2_gt_train'] <- summary(model4)$r.squared
      } else {
        df$prs_g <- prs_g[, i]
        df$prs_gtt <- prs_gtt[, i]
        model1 <- lm(df$res ~ prs_g[, i] + prs_gtt[, i])
        tunning[i, 'R2'] <- summary(model1)$r.squared
        model2 <- lm(df$res ~ prs_g[, i])
        tunning[i, 'R2_g'] <- summary(model2)$r.squared
        model3 <- lm(df$res ~ prs_gtt[, i])
        tunning[i, 'R2_gt'] <- summary(model3)$r.squared
        model_g <- lm(paste0('Y ~ prs_g + ', covar), data = df)
        df$res_g <- model_g$residuals
        model4 <- lm(res_g ~ prs_gtt, data = df)
        tunning[i, 'condR2_gt'] <- summary(model4)$r.squared

        df_train$prs_g <- prs_g_train[, i]
        df_train$prs_gtt <- prs_gtt_train[, i]
        model1 <- lm(df_train$res ~ prs_g_train[, i] + prs_gtt_train[, i])
        tunning[i, 'R2_train'] <- summary(model1)$r.squared
        model2 <- lm(df_train$res ~ prs_g_train[, i])
        tunning[i, 'R2_g_train'] <- summary(model2)$r.squared
        model3 <- lm(df_train$res ~ prs_gtt_train[, i])
        tunning[i, 'R2_gt_train'] <- summary(model3)$r.squared
        model_g <- lm(paste0('Y ~ prs_g + ', covar), data = df_train)
        df_train$res_g <- model_g$residuals
        model4 <- lm(res_g ~ prs_gtt, data = df_train)
        tunning[i, 'condR2_gt_train'] <- summary(model4)$r.squared
      }
    }
    tunning <- tunning[, c("lambda", "factor", "steps", "R2_g_train", "R2_gt_train", "condR2_gt_train", "R2_train", "R2_g", "R2_gt", "condR2_gt", "R2")]
    return(tunning)
  }

  R2_1 <- collect(1)
  R2_2 <- collect(2)
  R2_3 <- collect(3)
  R2_4 <- collect(4)

  R2 <- as.data.frame(cbind(R2_1[, 'R2'], R2_2[, 'R2'], R2_3[, 'R2'], R2_4[, 'R2']))
  mean_r2 <- apply(R2, 1, mean)

  if (fixG) {
    compare <- as.data.frame(cbind(R2_1[, 1:3], mean_r2))
    best <- compare[which.max(compare$mean_r2), ]
    rownames(best) <- NULL
    if (verbose) {
        cat("best_R2:", "\n")
        print(best)
    }
    s2_results <- list("best_R2" = best)
    #write.table(best, paste0(s2_output, '_initial_', initial, '_best_R2.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
  } else {
    condR2 <- as.data.frame(cbind(R2_1[, 'condR2_gt'], R2_2[, 'condR2_gt'], R2_3[, 'condR2_gt'], R2_4[, 'condR2_gt']))
    mean_cond_r2 <- apply(condR2, 1, mean)
    compare <- as.data.frame(cbind(R2_1[, 1:3], mean_r2, mean_cond_r2))
    
    best1 <- compare[which.max(compare$mean_r2), ]
    best2 <- compare[which.max(compare$mean_cond_r2), ]

    rownames(best1) <- NULL
    rownames(best2) <- NULL
    
    if (verbose) {
        cat("best_R2:", "\n")
        print(best1)
        cat("best_condR2:", "\n")
        print(best2)
    }
    
    s2_results = list("best_R2" = best1, "best_condR2" = best2)
    #write.table(best1, paste0(s2_output, '_initial_', initial, '_best_R2.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
    #write.table(best2, paste0(s2_output, '_initial_', initial, '_best_condR2.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
  }
  return(s2_results)
}

