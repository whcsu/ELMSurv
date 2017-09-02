## ---- eval=FALSE, message=FALSE, warning=FALSE, include=FALSE, results='asis'----
#  set.seed(123)
#  require(ELMSurv)
#  require(survival)
#  ## Survival Ensemble of ELM  with default settings
#  #Lung DATA
#  data(lung)
#  lung=na.omit(lung)
#  lung[,3]=lung[,3]-1
#  n=dim(lung)[1]
#  L=sample(1:n,ceiling(n*0.5))
#  trset<-lung[L,]
#  teset<-lung[-L,]
#  rii=c(2,3)
#  elmsurvfit=ELMSurvEN(x=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]),testx=teset[,-c(rii)])
#  # Get the 5th base model
#  basemodel=elmsurvfit[[1]]
#  #Print the c-index values
#  #library(survcomp)
#  #ci_elm=concordance.index(-rowMeans(elmsurvfit$precitedtime),teset$days,teset$status)[1]
#  #print(ci_elm)

