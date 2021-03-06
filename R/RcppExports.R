# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' A Buckley-James Esitmator
#' @title ELMSurv bjimpute
#' @param y  Survival time of training data. 
#' @param cen  The censoring indicator of training data.
#' @param x  The covariates(predictor variables) of training data.
#' @param inibeta   Initial values set for Buckley-James Imputation.
#' @return Imputed survival times
#' @seealso \code{\link{elm_surv}}
#' @author Hong Wang
#' @references
#' \itemize{
#'   \item Hong Wang et al (2017). A Survival Ensemble of Extreme Learning Machine. Applied Intelligence, in press.
#'  }
bjimpute <- function(y, cen, x, inibeta) {
    .Call(`_ELMSurv_bjimpute`, y, cen, x, inibeta)
}

#' A Kermatrix Compuation
#' @title ELMSurv kernmat
#' @param xtrain  The covariates(predictor variables) of training data.
#' @param kernel_type  Type of kernel matrix. kerneltype=1,a RBF kernel;kerneltype=2 , a linear kernel;kerneltype=3 ,a polynomial kernel;kerneltype=4, a sigmoid kernel.
#' @param kernel_para  Parameters for different types of kernels. A single value for kerneltype=1 or 2. A vector for kerneltype=3 or 4.
#' @param xtest   The covariates(predictor variables) of test data.
#' @return A kernel matrix
#' @seealso \code{\link{elm_surv}}
#' @author Hong Wang
#' @references
#' \itemize{
#'   \item Hong Wang et al (2017). A Survival Ensemble of Extreme Learning Machine. Applied Intelligence, in press.
#'  }
kernmat <- function(xtrain, kernel_type, kernel_para, xtest) {
    .Call(`_ELMSurv_kernmat`, xtrain, kernel_type, kernel_para, xtest)
}

