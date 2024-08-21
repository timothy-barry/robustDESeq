/*
#include <Rcpp.h>
#include <boost/random.hpp>
#include <boost/math/distributions/binomial.hpp>
using boost::math::binomial;
#include <stdexcept>
#include <algorithm>
#include "utilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
double compute_test_statistic(NumericVector a, NumericVector w, NumericMatrix D, const std::vector<int>& trt_idxs, int n_trt) {
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow();

  // iterate over the rows of D
  for (int i = 0; i < D_nrow; i ++) {
    inner_sum = 0;
    for (int j = 0; j < n_trt; j ++) {
      inner_sum += D(i, trt_idxs[j]);
    }
    lower_right += inner_sum * inner_sum;
  }

  // second, compute the lower-left hand of the denominator; also, compute the top
  for (int j = 0; j < n_trt; j ++) {
    top += a[trt_idxs[j]];
    lower_left += w[trt_idxs[j]];
  }

  // finally, compute the z-score
  return(top/sqrt(lower_left - lower_right));
}


// [[Rcpp::export]]
List run_adaptive_permutation_test(List precomp_list, IntegerVector x, int h, double alpha) {
  // define variables and objects
  int n = x.length(), m = precomp_list.length(), max_n_losses_active_set;
  std::vector<bool> active_set(m, true), futility_set(m, false), rejected_set(m, false);
  std::vector<double> stop_times(m), original_statistics(m), p_values(m);
  std::vector<int> n_losses(m, 0);
  double curr_test_stat, t = 0, h_doub = static_cast<double>(h), m_doub = static_cast<double>(m), threshold, n_in_active_set;
  List curr_precomp;

  // populate the trt_idxs vector
  std::vector<int> trt_idxs;
  for (int i = 0; i < x.length(); i++) if (x[i] == 1) trt_idxs.push_back(i);
  int n_trt = trt_idxs.size();
  if (n_trt == 0) throw std::invalid_argument("Zero treatment units.");

  // compute the original test statistics
  for (int i = 0; i < m; i++) {
    curr_precomp = precomp_list(i);
    original_statistics[i] = compute_test_statistic(curr_precomp(0), curr_precomp(1), curr_precomp(2), trt_idxs, n_trt);
  }

  // define objects related to random permutations
  std::vector<int> random_samp(n);
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);
  std::vector<double> i_doub_array(n_trt);
  for (int i = 0; i < n_trt; i ++) i_doub_array[i] = static_cast<double>(i);

  // iterate through time
  // for (int k = 0; k < 50000; k++) {
  for (int k = 0; k < 500; k++) {
    // increment time
    t++;
    // generate a random permutation
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n, n_trt);
    for (int j = 0; j < n; j ++) Rcout << random_samp[j] << " ";
    Rcout << "\n";

    // iterate over hypotheses
    for (int i = 0; i < m; i ++) {
      // if in the active set, compute the test statistic and update gamma
      if (active_set[i]) {
        // compute the test statistic
        curr_precomp = precomp_list(i);
        curr_test_stat = compute_test_statistic(curr_precomp(0), curr_precomp(1), curr_precomp(2), random_samp, n_trt);
        // determine whether we have a loss; if so, increment n_losses
        n_losses[i] += (curr_test_stat >= original_statistics[i] ? 1.0 : 0.0);
      }
    }

    // move elements from active set into futility set
    max_n_losses_active_set = 0;
    n_in_active_set = 0;
    for (int i = 0; i < m; i++) {
      if (active_set[i]) {
        if (n_losses[i] == h) { // hit loss limit
          active_set[i] = false;
          futility_set[i] = true;
          stop_times[i] = t;
        } else {
          n_in_active_set ++;
          if (n_losses[i] > max_n_losses_active_set) max_n_losses_active_set = n_losses[i];
        }
      }
    }

    // check for rejection of active set hypotheses
    threshold = t + h_doub - (h_doub * m_doub)/(n_in_active_set * alpha);
    if (static_cast<double>(max_n_losses_active_set) <= threshold) {
      for (int i = 0; i < m; i ++) {
        if (active_set[i]) {
          active_set[i] = false;
          rejected_set[i] = true;
          stop_times[i] = t;
        }
      }
      break;
    }
  }
  // compute the p-values
  for (int i = 0; i < m; i ++) {
    p_values[i] = h_doub/(stop_times[i] - static_cast<double>(n_losses[i]) + h_doub);
  }

  return List::create(Named("p_values") = p_values,
                      Named("rejected") = rejected_set);
}
*/
