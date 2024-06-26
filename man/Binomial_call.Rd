% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Binomialcall.R
\name{Binomial_call}
\alias{Binomial_call}
\title{Title}
\usage{
Binomial_call(n_steps, S, K, sigma, r, t)
}
\arguments{
\item{n_steps}{Amount of steps to simulate}

\item{S}{Stock price at start}

\item{K}{Strike price}

\item{sigma}{Volatility}

\item{r}{Interest rate}

\item{t}{Years until maturity of call}
}
\value{
Price of a  European call option calculated using the binomial tree
}
\description{
Calculates the Price of a European call option using a RCPP implementation
}
