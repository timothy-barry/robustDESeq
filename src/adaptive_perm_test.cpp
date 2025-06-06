#include <Rcpp.h>
#include <boost/random.hpp>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "utilities.h"
using namespace Rcpp;


// [[Rcpp::export]]
List run_adaptive_permutation_test_v2(List precomp_list, IntegerVector x, int side_code, int h, double alpha, int max_iterations, std::string test_stat_str, List custom_permutation_list) {
  // define variables and objects
  int n = x.length(), m = precomp_list.length(), n_trt;
  std::vector<bool> active_set(m, true), futility_set(m, false), rejected_set(m, false), curr_left_losses(m, false), curr_right_losses(m, false);
  std::vector<double> stop_times(m), original_statistics(m), p_values(m), n_right_losses(m, 0), n_left_losses(m, 0), n_losses(m, 0);
  std::vector<int> trt_idxs;
  double curr_test_stat, t = 0, h_doub = static_cast<double>(h), m_doub = static_cast<double>(m), threshold, n_in_active_set, max_n_losses_active_set;
  bool problem_round;

  // select the test statistic
  double (*funct)(List, const std::vector<int>&, int) = nullptr;
  if (test_stat_str == "compute_score_stat") {
    funct = compute_score_stat;
  } else if (test_stat_str == "compute_mean_over_treated_units") {
    funct = compute_mean_over_treated_units;
  } else if (test_stat_str == "compute_mw_test_statistic") {
    funct = compute_mw_test_statistic;
  } else {
    throw std::invalid_argument("Test statistic not recognized.");
  }

  // populate the trt_idxs vector
  for (int i = 0; i < x.length(); i++) if (x[i] == 1) trt_idxs.push_back(i);
  n_trt = trt_idxs.size();
  if (n_trt == 0) throw std::invalid_argument("Zero treatment units.");

  // compute the original test statistics
  for (int i = 0; i < m; i++) {
    original_statistics[i] = funct(precomp_list(i), trt_idxs, n_trt);
  }

  // define objects related to random permutations
  bool use_custom_permutations = custom_permutation_list.length() >= 1;
  boost::random::mt19937 generator(4);
  double n_doub;
  std::vector<int> random_samp;
  std::vector<double> i_doub_array;
  boost::random::uniform_real_distribution<double> distribution_real;
  boost::random::uniform_int_distribution<int> distribution_int;
  std::vector<std::vector<int>> custom_permutation_list_cpp;

  if (!use_custom_permutations) {
    distribution_real = boost::random::uniform_real_distribution<double>(0, 1);
    n_doub = static_cast<double>(n);
    i_doub_array = std::vector<double>(n_trt);
    for (int i = 0; i < n_trt; i++) i_doub_array[i] = static_cast<double>(i);
    random_samp = std::vector<int>(n);
    for (int i = 0; i < n; i++) random_samp[i] = i;
  } else {
    for (int i = 0; i < custom_permutation_list.length(); i++) {
      custom_permutation_list_cpp.push_back(Rcpp::as<std::vector<int>>(custom_permutation_list(i)));
    }
    int upper_value = custom_permutation_list.length() - 1;
    distribution_int = boost::random::uniform_int_distribution<int>(0, upper_value);
  }

  // iterate through time
  for (int k = 0; k < max_iterations; k++) {
    // set problem_round flag to false
    problem_round = false;
    // generate a random permutation
    if (!use_custom_permutations) {
      draw_wor_sample(generator, distribution_real, i_doub_array, random_samp, n_trt, n_doub);
    } else {
      random_samp = custom_permutation_list_cpp[distribution_int(generator)];
    }
    // iterate over hypotheses
    for (int i = 0; i < m; i ++) {
      // if in the active set, compute the test statistic
      if (active_set[i]) {
        // compute the test statistic
        curr_test_stat = funct(precomp_list(i), random_samp, n_trt);
        if (!std::isfinite(curr_test_stat)) {
          problem_round = true;
          break;
        }
        curr_left_losses[i] = (curr_test_stat <= original_statistics[i]);
        curr_right_losses[i] = (curr_test_stat >= original_statistics[i]);
      }
    }

    if (problem_round) {
      continue; // jump back to drawing wor sample
    } else {
      // increment time
      t++;
      // update n_losses (consider condensing this part)
      for (int i = 0; i < m; i++) {
        if (active_set[i]) {
          if (curr_right_losses[i]) n_right_losses[i] ++;
          if (curr_left_losses[i]) n_left_losses[i] ++;
        }
      }
      if (side_code == -1) {
        n_losses = n_left_losses;
      } else if (side_code == 0) {
        for (int i = 0; i < m; i ++) n_losses[i] = std::min(n_left_losses[i], n_right_losses[i]);
      } else {
        n_losses = n_right_losses;
      }
    }

    // move elements from active set into futility set
    max_n_losses_active_set = 0;
    n_in_active_set = 0;
    for (int i = 0; i < m; i++) {
      if (active_set[i]) {
        if (n_losses[i] >= h_doub) { // hit loss limit
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
    threshold = t + h_doub - ((side_code == 0 ? 2.0 : 1.0) * h_doub * m_doub)/(n_in_active_set * alpha);
    if (max_n_losses_active_set <= threshold) {
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
  // if any hypotheses remain in the active set, (implicitly) move to the futility set and update the stop time
  if (n_in_active_set >= 1) {
    for (int i = 0; i < m; i ++) {
      if (active_set[i]) {
        stop_times[i] = t;
      }
    }
  }

  // compute the p-values
  for (int i = 0; i < m; i ++) {
    p_values[i] = std::min(1.0, (side_code == 0 ? 2.0 : 1.0) * h_doub/(stop_times[i] - static_cast<double>(n_losses[i]) + h_doub));
  }

  return List::create(Named("p_values") = p_values,
                      Named("rejected") = rejected_set,
                      Named("n_losses") = n_losses,
                      Named("stop_times") = stop_times);
}
