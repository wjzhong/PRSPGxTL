
# function to adjust the phenotype 
adjust_y <- function(data, valid_gp, covar, fixG = FALSE){
    if (!is.null(valid_gp)) {
        data = data[data$group != valid_gp, ]
    }
    if (fixG) {
        model0=lm(paste0('Y ~ prs_g + ',covar), data=data)
    } else{
        model0=lm(paste0('Y ~ ',covar), data=data)
    }
    data$res = model0$residuals
    if (fixG) {
        model1=lm(res ~ prs_gt, data=data)
    } else {
        model1=lm(res ~ prs_g + prs_gt, data=data)
    }
    y_res = data[,c('ID','res','Tr')]
    if (fixG) {
        coeff = c(summary(model1)$coefficients['prs_gt','Estimate'])
    } else {
        coeff = c(summary(model1)$coefficients['prs_g','Estimate'], summary(model1)$coefficients['prs_gt','Estimate'])
    }
    return(list(y_res=y_res,coeff = coeff))
}
