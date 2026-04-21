#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Index_selection(const NumericMatrix& A, double b) {
  int n = A.nrow();
  List result(n);

  if (n == 0 || b <= 0) return result;
  int left = 0;
  int right = 0;
  for (int i = 0; i < n; ++i) {
    // Maintain left pointer: A[i,5] - A[left,5] < b
    while (left < n && A(i,4) - A(left,4) >= b) {
      left++;
    }
    // Maintain right pointer: A[right,5] - A[i,5] < b
    if (right < i) right = i;
    while (right < n && A(right,4) - A(i,4) < b) {
      right++;
    }
    // Count valid neighbors
    int count = 0;
    for (int j = left; j < right; ++j) {
      if (j == i) continue;
      bool cond12 =
        (A(i,0) != A(j,0)) ||
        (A(i,1) != A(j,1));
      bool cond34 =
        (A(i,2) != A(j,2)) ||
        (A(i,3) != A(j,3));
      if (cond12 && cond34) {
        count++;
      }
    }
    // Allocate and fill
    IntegerVector idx(count);
    int pos = 0;
    for (int j = left; j < right; ++j) {
      if (j == i) continue;
      bool cond12 =
        (A(i,0) != A(j,0)) ||
        (A(i,1) != A(j,1));
      bool cond34 =
        (A(i,2) != A(j,2)) ||
        (A(i,3) != A(j,3));
      if (cond12 && cond34) {
        idx[pos++] = j + 1; // R is 1-based
      }
    }
    result[i] = idx;
  }
  return result;
}
