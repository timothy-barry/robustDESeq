#include <boost/random.hpp>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

void draw_wor_sample(boost::random::mt19937& generator,
                     boost::random::uniform_real_distribution<double>& distribution,
                     const std::vector<double>& i_doub_array,
                     std::vector<int>& random_samp,
                     int n_trt, double n_doub) {
  double u;
  int pos;
  for (int i = 0; i < n_trt; i++) {
    u = distribution(generator);
    pos = floor((n_doub - i_doub_array[i]) * u);
    std::swap(random_samp[pos], random_samp[i]);
  }
  return;
}


// [[Rcpp::export]]
std::vector<std::vector<int>> generate_wor_sample_test(int n, int n_trt, int n_samples) {
  std::vector<std::vector<int>> out(n_samples);
  std::vector<int> random_samp(n);
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);
  double n_doub = static_cast<double>(n);
  std::vector<double> i_doub_array(n_trt);
  for (int i = 0; i < n_trt; i++) i_doub_array[i] = static_cast<double>(i);
  for (int i = 0; i < n; i++) random_samp[i] = i;
  for (int i = 0; i < n_samples; i ++) {
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n_trt, n_doub);
    out[i] = std::vector<int>(random_samp.begin(), random_samp.begin() + n_trt);
  }
  return out;
}
