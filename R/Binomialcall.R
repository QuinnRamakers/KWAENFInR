#' Title
#' @description
#'  Calculates the Price of a European call option using a RCPP implementation 
#' @param n_steps Amount of steps to simulate
#' @param S Stock price at start
#' @param K Strike price
#' @param sigma Volatility
#' @param r Interest rate
#' @param t Years until maturity of call
#'
#' @return Price of a  European call option calculated using the binomial tree
#' @export
#'
#' 
#' 
Binomial_call<-function(n_steps,S,K,sigma,r,t){
  if(t<=0){
    stop('Please provide a positive value for the time duration')
  }
  if(n_steps<=0){
    stop('Please provide a positive value for the amount of steps')
  }
 return( .B_numeraire_Call(n_steps,S,K,sigma,r,t))
}