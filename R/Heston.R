
#' Title
#' @description
#' Function that simulates a Heston model using Monte-Carlo
#' 
#' @param rho Correlation for the Brownian motions
#' @param sigma Variance parameter for the Vt process
#' @param lambda Mean-reversion rate parameter
#' @param t Amount of periods to simulate
#' @param r Interest rate
#' @param theta Long term mean for Vt
#' @param V0 Value of Vt at the start
#' @param x0 Initial log of stock price
#' @param n_steps Amount of time steps 
#' @param paths Amount of paths to simulate 
#'
#' @return HestonMc object where every row is a simulation and every column a timestep
#' @export
#'
#' @examples
HESTON<-function(rho,sigma,lambda,t,r,theta,V0,x0,n_steps,paths){
  dt=t/n_steps
  x_price=matrix(0,nrow=paths,ncol=n_steps+1)
  stocvol=matrix(0,nrow=paths,ncol=n_steps+1)
  x_price[,1]=x0
  stocvol[,1]=V0
  Sigma_mat <- matrix(c(1, rho,
                        rho, 1), 
                      nrow = 2, ncol = 2, byrow = TRUE)
  for(i in 2:(n_steps+1)){
    Wt=mvtnorm::rmvnorm(paths, sigma = Sigma_mat)
    x_price[,i]=x_price[,i-1]+(r-0.5*stocvol[,i-1])*dt+sqrt(stocvol[,i-1]*dt)*Wt[,1]
    stocvol[,i]=stocvol[,i-1]+lambda*(theta-stocvol[,i-1])*dt+sigma*sqrt(stocvol[,i-1]*dt)*Wt[,2]
    #set negative volatility to 0
    stocvol[stocvol[,i]<0,i]=abs(stocvol[stocvol[,i]<0,i])
  }
  systemout=HESTON_mc(x_price,stocvol)
  
  return(systemout)
}

HESTON_mc <- function(prices,volatility) {
  obj <- list(prices = prices,volatility=volatility )
  class(obj) <- "Heston_mc"
  return(obj)
}
plot.Heston_mc<-function(x, ...){
  mean_prices_df <- data.frame(Time = 0:(ncol(x$prices) - 1), value = colMeans(x$prices), Path = "Mean")
  mean_volatility_df <- data.frame(Time = 0:(ncol(x$volatility) - 1), value = colMeans(x$volatility), Path = "Mean")
  mean_prices_df$Type <- "Price"
  mean_volatility_df$Type <- "Volatility"
  combined_df <- rbind(mean_prices_df, mean_volatility_df)
  return(ggplot(combined_df, aes(x = Time, y = value, group = Path, color = Type)) +
    geom_line(size = 1.5) +
    labs(
      title = "Development of Mean Price and Volatility",
      subtitle = "Simulated Mean Price and Volatility Paths Over Time",
      x = "Timestep",
      y = "Value"
    ) +
    facet_wrap(vars(Type),scales = 'free_y')
  +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    ))
}
