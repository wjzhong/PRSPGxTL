
#' @title inner_layer_CV
#' @description Inner layer of cross-validation.
#' @param sum_stats dataframe including PRS weights from disease GWAS.
#' @param ped dataframe including phenotype and covariates.
#' @param G The genotype matrix of individual x SNP. Row names are individual ID and column names are SNP ID. 
#' @param bim Allele info for the SNPs in genotype file. The following columns are required: SNP, Ref. Column names are needed.
#' @param initial Indicating the starting values of the predictive effects, 'PRS' or 'zero'.
#' @param num_snp number of total snps with non-zero effects. 
#' It is used to determine the learning rate as 
#' (factor * 1 / num_snp), where the factor is varied across a grid of values: 1, 10, 50, 70. 
#' @param covar covariates to be adjusted in the model (other than Tr), such as "age,gender". This is optional. Default is NULL. 
#' @param fixG logical, indicating whether the genetic main (G) effect is fixed. Default is FALSE.
#' @param verbose logical, indicating whether to print more details. Default is FALSE.
#' @return A list of S_prog and S_pred values for the validation set, the training set, and the tuning parameters. 
#' @export
#' @examples
#' \dontrun{
#' data(sum_stats)
#' data(ped)
#' data(G)
#' data(bim)
#' result_list = inner_layer_CV(sum_stats = sum_stats, ped = ped, G = G, bim = bim, 
#'  initial = "PRS", num_snp = 2349, covar = NULL, fixG = FALSE, verbose = TRUE)
#' }

inner_layer_CV <- function(sum_stats, ped, G, bim, initial, num_snp, covar = NULL, fixG = FALSE, verbose = FALSE) {

    if (length(covar) == 0) {
        covar <- 'Tr'
    } else {
        covars <- unlist(strsplit(covar, ","))
        covar <- paste0(c('Tr', covars), collapse = " + ")
    }

    ped$prs_gt = ped$Tr * ped$prs_g
    G = G[rownames(G) %in% ped$ID, ,drop=FALSE]

    ######### 0. match snp and flip +/-
    bim2 = bim[,c('SNP','Ref')]
    colnames(bim2)[2] = 'Ref_G'

    bim_sum_stats=inner_join(sum_stats, bim2, by="SNP")
    bim_sum_stats$Beta2=NA

    #ref in prs matched with ref in training
    flag1=which(bim_sum_stats$Ref==bim_sum_stats$Ref_G)
    if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
    #eff in prs matched with ref in training
    flag2=which(bim_sum_stats$Eff==bim_sum_stats$Ref_G)
    if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2]}
    #table(bim_sum_stats$Beta2 == bim_sum_stats$Beta_matched)
    bim_sum_stats=bim_sum_stats[is.na(bim_sum_stats$Beta2)==F,]

    
    ######### split into 4 (3 for training, 1 for validation)
    group_assignment <- rep(1:4, each = nrow(ped) %/% 4)
    if (nrow(ped) %% 4 > 0){ group_assignment <- c(group_assignment, 1:(nrow(ped) %% 4))}
    set.seed(123)  
    # Randomly assign each individual to one of 4 groups
    ped$group <- sample(group_assignment)
    group_values <- as.numeric(names(table(ped$group)))
    #valid_gp = group_values[1]

    result_list = list()
    for (j in 1:length(group_values)){
        if (verbose) {
            print(paste0("Inner CV: ", j))
        }
        valid_gp = group_values[j]
        # adjust the phenotype 
        obj = adjust_y(ped, valid_gp, covar, fixG)
        # 
        if (fixG) {
            step = 200
        } else {
            step = 30
        }
        # Gradient descent and keep all beta
        GD_obj = GD(obj, G, bim_sum_stats, step, initial, num_snp, NULL, NULL, fixG)
        # apply all the beta on the validation to get the scoresum (for this region)
        result_list[[j]] = get_score(GD_obj, obj, G, fixG)
        #save(result, file=paste0(output,'_initial_',initial,'_innerCV_',j,'.Rdata'))
    }
    return(result_list)
}


# function to apply all the beta on the validation to get the scoresum (for this region)
get_score <- function(GD_obj, obj, G, fixG = FALSE) {
    oc = GD_obj$oc
    y_res = obj$y_res
    
    if (fixG) {
        oc_gt = oc
        snps = sub("_gt$", "", rownames(oc_gt))
    } else {
        oc_g = oc[1:(nrow(oc) / 2), ]
        oc_gt = oc[(nrow(oc) / 2 + 1):nrow(oc), ]
        snps = rownames(oc_g)
    }

    geno_valid = G[!rownames(G) %in% y_res$ID, snps]
    geno_train = G[rownames(G) %in% y_res$ID, snps]

    # Check consistency
    if (!fixG) {
        stopifnot(all(colnames(geno_valid) == rownames(oc_g)))
    }

    # Convert to matrix
    geno_valid = as.matrix(geno_valid)
    geno_train = as.matrix(geno_train)
    
    if (fixG) {
        weight_gt = oc_gt[, 6] / oc_gt$sd
        score_gt = geno_valid %*% weight_gt
        score_gt_train = geno_train %*% weight_gt

        for (j in 7:ncol(oc_gt)) {
            weight_gt = oc_gt[, j] / oc_gt$sd
            score_gt = cbind(score_gt, geno_valid %*% weight_gt)
            score_gt_train = cbind(score_gt_train, geno_train %*% weight_gt)
        }
        
        return(list(score_gt = score_gt, score_gt_train = score_gt_train, tunning = GD_obj$tunning))
    } else {
        weight_g = oc_g[, 6] / oc_g$sd
        weight_gt = oc_gt[, 6] / oc_gt$sd
        score_g = geno_valid %*% weight_g
        score_gt = geno_valid %*% weight_gt
        score_g_train = geno_train %*% weight_g
        score_gt_train = geno_train %*% weight_gt

        for (j in 7:ncol(oc_g)) {
            weight_g = oc_g[, j] / oc_g$sd
            weight_gt = oc_gt[, j] / oc_gt$sd
            score_g = cbind(score_g, geno_valid %*% weight_g)
            score_gt = cbind(score_gt, geno_valid %*% weight_gt)
            score_g_train = cbind(score_g_train, geno_train %*% weight_g)
            score_gt_train = cbind(score_gt_train, geno_train %*% weight_gt)
        }
        
        return(list(score_g = score_g, score_gt = score_gt, score_g_train = score_g_train, score_gt_train = score_gt_train, tunning = GD_obj$tunning))
    }
}


