#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Recursive backward induction function
arma::vec B_backward_induction(const arma::vec& current_layer, double p_star) {
  if (current_layer.n_elem == 1) {
    return current_layer * 100;
  } else {
    arma::vec previous_layer = p_star * current_layer.subvec(0, current_layer.n_elem - 2) +
      (1 - p_star) * current_layer.subvec(1, current_layer.n_elem - 1);
    return B_backward_induction(previous_layer, p_star);
  }
}

// [[Rcpp::export(.B_numeraire_Call)]]
double B_numeraire_Call(int n_steps, double S, double K, double sigma, double r, double T) {
  double delta = T / n_steps;
  double R = exp(r * delta);
  double d = R * exp(-sigma * sqrt(delta));
  double u = R * exp(sigma * sqrt(delta));
  double p_star = (R * S - d * S) / (u * S - d * S);
  
  arma::vec domain = arma::regspace<arma::vec>(0, n_steps);
  arma::vec d_powers = arma::pow(arma::ones<arma::vec>(domain.n_elem) * d, domain);
  arma::vec u_powers = arma::pow(arma::ones<arma::vec>(domain.n_elem) * u, n_steps - domain);
  arma::vec maturity_value = arma::max(S * d_powers % u_powers, K * arma::ones<arma::vec>(domain.n_elem)) / (S * exp(r * T));
  
  arma::vec numeraire = S * arma::exp(domain * log(R));
  
  arma::vec result = B_backward_induction(maturity_value, p_star);
  return result(0); // Return the first element as the final result
}
