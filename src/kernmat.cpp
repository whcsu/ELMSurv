#include "global.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//' A Kermatrix Compuation
//' @title ELMSurv kernmat
//' @param xtrain  The covariates(predictor variables) of training data.
//' @param kernel_type  Type of kernel matrix. kerneltype=1,a RBF kernel;kerneltype=2 , a linear kernel;kerneltype=3 ,a polynomial kernel;kerneltype=4, a sigmoid kernel.
//' @param kernel_para  Parameters for different types of kernels. A single value for kerneltype=1 or 2. A vector for kerneltype=3 or 4.
//' @param xtest   The covariates(predictor variables) of test data.
//' @return A kernel matrix
//' @seealso \code{\link{elm_surv}}
//' @author Hong Wang
//' @references
//' \itemize{
//'   \item Hong Wang et al (2017). A Survival Ensemble of Extreme Learning Machine. Applied Intelligence, in press.
//'  }
//[[Rcpp::export]]
Rcpp::NumericMatrix  kernmat( SEXP xtrain,SEXP kernel_type, SEXP kernel_para,SEXP xtest) {
 //kernel_type 1:RBF_kernel 2:lin_kernel 3:poly_kernel 4:sigmoid_kernel 
 int kerneltype=as <int> (kernel_type);
 Rcpp::NumericVector  kernel_pa(kernel_para);
 
 //convert a dataframe to matrix
	Rcpp:: NumericMatrix xtrNM=testDFtoNM(xtrain);
	int ntr = xtrNM.nrow();
	
    if (kerneltype==1){
      if (Rf_isNull(xtest)){
		Rcpp:: NumericMatrix xhh(ntr,ntr);
        for (int i=0;i<ntr;i++)
          	for (int j=0;j<ntr;j++)
			{
			 Rcpp::NumericVector diff = xtrNM(i,_)-xtrNM(j,_);
			 Rcpp::NumericVector diff2=pow(diff,2);
			 //std::cout<<diff[0];
			  xhh(i,j)=std::accumulate(diff2.begin(), diff2.end(), 0.0);		  
				  
			} 
		
		   for (int i=0;i<ntr;i++){
			      //exp support only vectors not matrix
				  xhh(i,_) =exp(-xhh(i,_)/kernel_pa[0]);
				   //std::cout<<xhh(i,0)<<"\n";		
		   }
		return(xhh);
      }
       else{
		 Rcpp:: NumericMatrix xteNM=testDFtoNM(xtest);
	     int nte = xteNM.nrow();               
         Rcpp:: NumericMatrix xhh(nte,ntr);
		
		 
         for (int i=0;i<nte;i++)
          	for (int j=0;j<ntr;j++)
			{
			  Rcpp::NumericVector diff = xteNM(i,_)-xtrNM(j,_);
			   Rcpp::NumericVector diff2=pow(diff,2);
			   xhh(i,j)=std::accumulate(diff2.begin(), diff2.end(), 0.0);	
			} 
		
		   for (int i=0;i<nte;i++){
				  xhh(i,_) =exp(-xhh(i,_)/kernel_pa[0]);
				  //std::cout<<xhh(i,0)<<"\n";	
		   }
       return(xhh);
      }
    }

	if (kerneltype==2){
      if (Rf_isNull(xtest)){
		  //wrap:: From arma::mat to Rcpp::NumericMatrix
		Rcpp:: NumericMatrix xhh=wrap(mm_transpose(xtrNM,xtrNM)); 
         for (int i=0;i<ntr;i++){
				  xhh(i,_) =xhh(i,_)+kernel_pa[0];				 
		   }
		
		return(xhh);
      }
       else{
		 Rcpp:: NumericMatrix xteNM=testDFtoNM(xtest);
	     int nte = xteNM.nrow();               
         Rcpp:: NumericMatrix xhh=wrap(mm_transpose(xteNM,xtrNM)); 
		
		  for (int i=0;i<nte;i++){
				  xhh(i,_) =xhh(i,_)+kernel_pa[0];	
				  //std::cout<<xhh(i,0)<<"\n";	
		   }
       return(xhh);
      }
    }
	
	if (kerneltype==3 && kernel_pa.size()>1){

      if (Rf_isNull(xtest)){
		  //wrap:: From arma::mat to Rcpp::NumericMatrix
		Rcpp:: NumericMatrix xhh=wrap(mm_transpose(xtrNM,xtrNM)); 
         for (int i=0;i<ntr;i++){
				  xhh(i,_)=xhh(i,_)+kernel_pa[0];
                  xhh(i,_) =pow(xhh(i,_),kernel_pa[1]);				  
		   }
		
		return(xhh);
      }
       else{
		 Rcpp:: NumericMatrix xteNM=testDFtoNM(xtest);
	     int nte = xteNM.nrow();               
         Rcpp:: NumericMatrix xhh=wrap(mm_transpose(xteNM,xtrNM)); 
		
		  for (int i=0;i<nte;i++){
				  xhh(i,_)=xhh(i,_)+kernel_pa[0];
                  xhh(i,_) =pow(xhh(i,_),kernel_pa[1]);		
				  //std::cout<<xhh(i,0)<<"\n";	
		   }
       return(xhh);
      }
    }
	
	if (kerneltype==4 && kernel_pa.size()>1){

      if (Rf_isNull(xtest)){
		  //wrap:: From arma::mat to Rcpp::NumericMatrix
		Rcpp:: NumericMatrix xhh=wrap(mm_transpose(xtrNM,xtrNM)); 
         for (int i=0;i<ntr;i++){
			xhh(i,_)=xhh(i,_)*kernel_pa[0]+kernel_pa[1];
            xhh(i,_) =(exp(xhh(i,_))-exp(-xhh(i,_)))/(exp(xhh(i,_))+exp(-xhh(i,_)));				  
		   }
		
		return(xhh);
      }
       else{
		 Rcpp:: NumericMatrix xteNM=testDFtoNM(xtest);
	     int nte = xteNM.nrow();               
         Rcpp:: NumericMatrix xhh=wrap(mm_transpose(xteNM,xtrNM)); 
		
		  for (int i=0;i<nte;i++){
			xhh(i,_)=xhh(i,_)*kernel_pa[0]+kernel_pa[1];
            xhh(i,_) =(exp(xhh(i,_))-exp(-xhh(i,_)))/(exp(xhh(i,_))+exp(-xhh(i,_)));	
		   }
       return(xhh);
      }
    }
	
	
return R_NilValue;		
	
}