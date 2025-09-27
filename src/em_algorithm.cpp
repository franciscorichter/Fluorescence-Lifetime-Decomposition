#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List EMAlgorithmCpp(IntegerVector pixel_indices, NumericVector tau_vals,
                    NumericMatrix pi_hat, NumericVector lambda_hat,
                    IntegerMatrix coords, int n_rows, int n_cols, int max_iter) {
  int N = tau_vals.size();
  int N_pix = pi_hat.nrow();
  int m = pi_hat.ncol();

  NumericVector loglikelihoods(max_iter);
  NumericMatrix lambda_hist(max_iter, m);

  NumericMatrix resp(N, m);
  NumericMatrix numerators(N, m);
  NumericVector row_sum(N);

  for (int iter = 0; iter < max_iter; iter++) {
    // E-step: compute responsibilities
    for (int i = 0; i < N; i++) {
      row_sum[i] = 0;
      for (int d = 0; d < m; d++) {
        int pix = pixel_indices[i] - 1; // convert to 0-index
        numerators(i, d) = pi_hat(pix, d) * lambda_hat[d] * exp(-lambda_hat[d] * tau_vals[i]);
        row_sum[i] += numerators(i, d);
      }
    }
    for (int i = 0; i < N; i++) {
      for (int d = 0; d < m; d++) {
        resp(i, d) = numerators(i, d) / row_sum[i];
      }
    }

    // Compute log-likelihood
    double logL = 0.0;
    for (int i = 0; i < N; i++) {
      logL += log(row_sum[i]);
    }
    loglikelihoods[iter] = logL;

    // M-step: update global decay rates
    for (int d = 0; d < m; d++) {
      double num = 0.0, den = 0.0;
      for (int i = 0; i < N; i++) {
        num += resp(i, d);
        den += resp(i, d) * tau_vals[i];
      }
      lambda_hat[d] = num / den;
    }

    // M-step: update pixel-specific mixing proportions using 3x3 neighborhood
    for (int pix = 0; pix < N_pix; pix++) {
      int i0 = coords(pix, 0);  // already 1-indexed
      int j0 = coords(pix, 1);
      int i_min = std::max(1, i0 - 1);
      int i_max = std::min(n_rows, i0 + 1);
      int j_min = std::max(1, j0 - 1);
      int j_max = std::min(n_cols, j0 + 1);
      NumericVector sum_resp(m);
      int count = 0;
      for (int p = 0; p < N_pix; p++) {
        int ip = coords(p, 0);
        int jp = coords(p, 1);
        if (ip >= i_min && ip <= i_max && jp >= j_min && jp <= j_max) {
          for (int i = 0; i < N; i++) {
            if (pixel_indices[i] - 1 == p) {
              for (int d = 0; d < m; d++) {
                sum_resp[d] += resp(i, d);
              }
              count++;
            }
          }
        }
      }
      if (count > 0) {
        for (int d = 0; d < m; d++) {
          pi_hat(pix, d) = sum_resp[d] / count;
        }
        double total = 0.0;
        for (int d = 0; d < m; d++) total += pi_hat(pix, d);
        for (int d = 0; d < m; d++) {
          pi_hat(pix, d) = std::max(pi_hat(pix, d), 1e-6) / total;
        }
      }
    }

    for (int d = 0; d < m; d++) {
      lambda_hist(iter, d) = lambda_hat[d];
    }

    Rcout << "Iteration " << iter+1 << ": Log-Likelihood = " << std::round(logL*100)/100.0
          << " | Rates = ";
    for (int d = 0; d < m; d++) {
      Rcout << std::round(lambda_hat[d]*1000)/1000.0 << (d < m-1 ? ", " : "\n");
    }
  }

  return List::create(Named("pi_hat") = pi_hat,
                      Named("lambda_hat") = lambda_hat,
                      Named("loglikelihoods") = loglikelihoods,
                      Named("lambda_hist") = lambda_hist);
}
