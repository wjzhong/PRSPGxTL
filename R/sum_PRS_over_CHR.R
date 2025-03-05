
#' @title sum_PRS_over_CHR
#' @description sum the PRS over all chromosomes for the testing samples 
#' @param s3_byCHR_files a vector or list of file names for output from "rerun_algorithm" function for each chromosome
#' @param fixG logical, indicating whether the genetic main (G) effect is fixed. Default is FALSE.
#' @return a data frame with the sum of PRS over all chromosomes for the testing samples 
sum_PRS_over_CHR <- function(s3_byCHR_files, fixG = FALSE) {
  if (fixG) {
    test <- readRDS(s3_byCHR_files[[1]])[["testScoresum"]]
    prs_gt <- test$score_gt
    
    for (i in 2:22) {
      test <- readRDS(s3_byCHR_files[[i]])[["testScoresum"]]
      prs_gt <- prs_gt + test$score_gt
    }
    
    score <- as.data.frame(prs_gt)
    colnames(score) <- c('prs_gt')
    score$ID <- rownames(score)
    score[, 1] <- scale(score[, 1])
    
  } else {
    
    test <- readRDS(s3_byCHR_files[[1]])[["testScoresum"]]
    prs_g <- test$score_g
    prs_gt <- test$score_gt
    
    for (i in 2:22) {
      test <- readRDS(s3_byCHR_files[[i]])[["testScoresum"]]
      prs_g <- prs_g + test$score_g
      prs_gt <- prs_gt + test$score_gt
    }
    
    score <- as.data.frame(cbind(prs_g, prs_gt))
    colnames(score)[1:2] <- c('prs_g', 'prs_gt')
    score$ID <- rownames(score)
    score[, 1:2] <- scale(score[, 1:2])
   
  }
  return(score)
}


