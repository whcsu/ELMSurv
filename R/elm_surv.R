##' Extreme Learning Machine Using the Buckley-James estimator
##' @title ELMSurv elm_surv
##' @param trainx  The covariates(predictor variables) of training data.
##' @param trainy  Survival time and censored status of training data. Must be a Surv  \code{survival} object
##' @param testx  The covariates(predictor variables) of test data.
##' @param Regularization_coefficient  Ridge or Tikhonov regularization parameter. Default value for \code{\link{ELMSurvEN}} is 10000. It need be set by the user here when using a single base ELM survival model. Also known as \eqn{C} in the ELM paper.
##' @param kerneltype Type of kernel matrix. kerneltype=1,a RBF kernel;kerneltype=2 , a linear kernel;kerneltype=3 ,a polynomial kernel;kerneltype=4, a sigmoid kernel.
##' @param Kernel_para Parameters for different types of kernels. A single value for kerneltype=1 or 2. A vector for kerneltype=3 or 4.
##' @return List of returned values
##'   \tabular{ll}{
##'       \code{elmsurvfit}    \tab  Mean Square Error(MSE) on training data. \cr
##'       \code{newy} \tab Esitmated survival times of training data by the Buckley-James estimator. \cr
##'       \code{outputWeight} \tab Weights of the output layer in ELM. \cr
##'       \code{testpre} \tab The estimated survival times for \code{testx} data. \cr
##'   }
##' @seealso \code{\link{ELMSurvEN}}
##' @author Hong Wang
##' @references
##' \itemize{
##'   \item Hong Wang et al (2017). A Survival Ensemble of Extreme Learning Machine. Applied Intelligence, in press.
##'  }
##' @export
elm_surv <-
function(trainx,trainy,testx,Regularization_coefficient, kerneltype, Kernel_para)
{ 

  ny <- ncol(trainy)

  status <- trainy[, ny]
  survtime = trainy[, 1L]
  
  #imputey
  
  newy = bjimpute(y = survtime, cen = status, x = trainx,inibeta = NULL)
  
  c = Regularization_coefficient
 
  omega_train = kernmat(trainx,kerneltype, Kernel_para,NULL)
  #Calculate the mode coefficient beta
  outputWeight = solve(omega_train + diag(rep(1, nrow(trainx)))/c) %*% newy
  #Calculate the training output
  trainypre = omega_train %*% outputWeight 
  
  trainMSE = sqrt(sum((trainypre - survtime)^2))
   
  #Calculate the output of testing input   
  #omega_test=kernel_matrix(trainx,Kernel_type,Kernel_para,testx)
  omega_test = kernmat(trainx,kerneltype, Kernel_para,testx)
  testpre = omega_test %*% outputWeight
  
  
  return(list(trainMSE = trainMSE,newy =newy, outputWeight = outputWeight,testpre = testpre))  
}
