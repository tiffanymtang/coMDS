// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;


// Get K nearest neighbors given distance vector x
// [[Rcpp::export]]
IntegerVector get_neighbors_vec(NumericVector x, double epsilon) {
  // Create a vector of pairs (value, index)
  int n = x.size();
  std::vector<std::pair<double, int>> x_idxs(n);
  for (int i = 0; i < n; ++i) {
    x_idxs[i] = std::make_pair(x[i], i); // Store value and its original index
  }

  // Use partial_sort to sort the pairs by value in ascending order
  std::sort(x_idxs.begin(), x_idxs.end(),
                    [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                      return a.first < b.first; // Sort in ascending order
                    });

  // Prepare output vectors for largest K values and their indices
  IntegerVector neighbor_idxs;
  for (int k = 0; k < n; ++k) {
    if(x_idxs[k].first < epsilon){
      neighbor_idxs.push_back(x_idxs[k].second);
    }
  }
  // Sort the indices in ascending order
  std::sort(neighbor_idxs.begin(), neighbor_idxs.end());

  return neighbor_idxs;
}

// Get pi-th percentile neighbors given a list of distance matrices (each encoded as a vector)
// [[Rcpp::export]]
List get_neighbors(List Ds, NumericVector epsilon_vec) {
  NumericVector temp = Ds[0];
  int n = (1 + sqrt(1 + 8 * temp.size())) / 2;
  int L = Ds.size();

  List neighbor_idxs_ls(L);
  for (int l = 0; l < L; ++l) {
    NumericVector Dl = Ds[l];
    List neighbor_idxs(n);
    double epsilon = epsilon_vec[l];
    for (int i = 0; i < n; ++i) {
      NumericVector Di(n);
      for (int j = 0; j < n; ++j) {
        if (i == j) {
          Di[j] = max(Dl) + 1;
        } else if (i < j) {
          Di[j] = Dl(i * (2 * n - i - 1) / 2 + (j - (i + 1)));
        } else {
          Di[j] = Dl(j * (2 * n - j - 1) / 2 + (i - (j + 1)));
        }
      }
      neighbor_idxs[i] = get_neighbors_vec(Di, epsilon);
    }
    neighbor_idxs_ls[l] = neighbor_idxs;
  }
  return neighbor_idxs_ls;
}
