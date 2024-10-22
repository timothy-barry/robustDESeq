#include <boost/random.hpp>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector compute_score_stat_benchmark(const List& precomp, int s, const List& trt_idxs_list) {
  NumericVector a = precomp(0);
  NumericVector w = precomp(1);
  List D_list = precomp(2);

  double lower_right, lower_left, top, inner_sum;
  int D_nrow = D_list.length();
  NumericVector d;
  int B = trt_idxs_list.length();
  NumericVector out(B);
  IntegerVector trt_idxs;

  for (int k = 0; k < B; k++) {
    trt_idxs = trt_idxs_list(k);
    lower_right = 0; lower_left = 0; top = 0;
    // iterate over the rows of D
    for (int i = 0; i < D_nrow; i ++) {
      inner_sum = 0;
      d = D_list(i);
      for (int j = 0; j < s; j ++) {
        inner_sum += d[trt_idxs[j]];
      }
      lower_right += inner_sum * inner_sum;
    }

    // second, compute the lower-left hand of the denominator; also, compute the top
    for (int j = 0; j < s; j ++) {
      top += a[trt_idxs[j]];
      lower_left += w[trt_idxs[j]];
    }

    // finally, compute the z-score
    out[k] = top/sqrt(lower_left - lower_right);
  }
  return out;
}
