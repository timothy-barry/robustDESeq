#include <Rcpp.h>
#include <boost/random.hpp>
using namespace Rcpp;

#ifndef DRAW_WOR_SAMPLE
#define DRAW_WOR_SAMPLE
void draw_wor_sample(boost::random::mt19937& generator, boost::random::uniform_real_distribution<double>& distribution, const std::vector<double>& i_doub_array, std::vector<int>& x, int n_tot, int M);
#endif
