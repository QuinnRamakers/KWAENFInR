
update_DF<-function(ft,dT,N,DFt,t){
  
  for( i in (t+1):(N+1) ){ DFt[i] <- DFt[i-1]/(1+dT*ft[i-1])}
  return(DFt)
}
#drift term calculation
drift_calc<-function(ft,dT,N,tdrift,sigma,t){
  
  #for the final year the for loop bounds switch   so this is needed to catch it
  if(t+1==N){
    return(-tdrift)
  }
  for( i in (N-1):(t+1) ){ tdrift[i] <- (tdrift[i+1]+(dT*sigma[i])/(1+dT*ft[i+1]))}
  return(-tdrift)
}
#method that pushes forward using euler discretisation 
dfi_calc<-function(ft,dT,N,tdrift,sigma,t){
  (ft+tdrift*sigma+sigma*rnorm(1))
}

#' Title
#' @description
#' This function simulates LIBOR market rate models and returns an object with 
#' the matrixes for multiple simulations to be used for pricing of desired derivative.
#' @details
#' We implement the LIBOR market rate model based on martingale pricing and compensate 
#' using the Radon-Nikodyn derivative to add a drift terms as found in the drift matrix. 
#' This method always simulates untill the highest maturity is reached.
#' 
#' @param forwardrates Vector of forward rates with the same frequency. 
#'
#' @param dT Timesteps where 1 is equal to a year. This should be equal to the time between forward rates. 
#' @param sigma Volatility matrix, in case of constant volatility provide a matrix with only a singular value for all indexes. 
#'
#' @return Object of class LIBORmc containing all matrixes simulated. 
#' @export
#'
#' @examples
LIBORSIM<- function(forwardrates,dT,sigma){
  #set matrices
  N=length(forwardrates)
  LIBOR=matrix(NA,N,N)
  drift=matrix(NA,N,N)
  DF=matrix(NA,N+1,N)
  #Starting State
  LIBOR[,1]=forwardrates
  #diagonal of df=1
  diag(DF)<-1
  #drift of TN=0 as it is the measure we use.
  drift[N,c(2:N)]=0
  
  DF[,1]=update_DF(LIBOR[,1],dT,N,DF[,1],1)
  for(i in 2:N){
    drift[,i]=drift_calc(LIBOR[,i-1],dT,N,drift[,i],sigma[,i-1],i-1)
    LIBOR[,i]=dfi_calc(LIBOR[,i-1],dT,N,drift[,i],sigma[,i-1],i-1)
    DF[,i]=update_DF(LIBOR[,i],dT,N,DF[,i],i)
  }
  return(list("LIBOR" = LIBOR, "DF" = DF,"Drift"=drift))
}


mcLIBOR<-function(forwardrates,dT,sigma,n){

 
  cl<-parallel::makeCluster(parallel::detectCores()-1)
  registerDoParallel(cl)
  mc<-foreach::foreach(b = 1:n,.combine='rbind',.export = c('LIBORSIM','dfi_calc','drift_calc','update_DF')) %dopar%{
    LIBORSIM(forwardrates,dT,sigma)
  }
  parallel::stopCluster(cl)
  result <- LIBORmc(rates = mc[,1], DF = mc[,2], drift = mc[,3])
  return(result)
}
LIBORmc <- function(rates,DF,drift) {
  obj <- list(rates = rates, DF = DF,drift = drift )
  class(obj) <- "LIBORmc"
  return(obj)
}
print.LIBORmc <- function(x, ...) {
  cat("LIBOR Monte Carlo Simulation Results\n")
  cat("Drift:\n")
  print(x$drift)
  cat("\nDiscount Factors (DF):\n")
  print(x$DF)
  cat("\nRates:\n")
  print(x$rates)
}
