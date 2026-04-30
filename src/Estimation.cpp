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

// Scalar version of k(r)
inline double k_cpp(double r) {
  if (std::abs(r) < 1.0) {
    return 0.75 * (1.0 - r * r);
  } else {
    return 0.0;
  }
}

// Scalar version of k_b(r, b)
// Useful inside tight C++ loops
inline double k_b_scalar(double r, double b) {
  if (b == 0.0) {
    return R_PosInf;
  } else if (r < b) {
    return k_cpp(r / b) / b;
  } else {
    return 0.0;
  }
}

// Vectorized version of k_b(r, b)
NumericVector k_b_cpp(NumericVector r, double b) {
  int n = r.size();
  NumericVector out(n);

  if (b == 0.0) {
    std::fill(out.begin(), out.end(), R_PosInf);
    return out;
  }

  for (int i = 0; i < n; ++i) {
    if (r[i] < b) {
      out[i] = k_cpp(r[i] / b) / b;
    } else {
      out[i] = 0.0;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector compute_c0_cpp(NumericVector dist,
                             NumericVector Z_v,
                             NumericVector e,
                             double lambda,
                             double b,
                             double N_tau) {

  int n = dist.size();
  NumericVector c0(n, NA_REAL);

  const double two_pi = 2.0 * M_PI;

  for (int i = 0; i < n; ++i) {

    int left = i;
    int right = i;

    // Move left boundary
    while (left > 0 && dist[i] - dist[left] < b) {
      --left;
    }

    // Move right boundary
    while (right < n - 1 && dist[right] - dist[i] < b) {
      ++right;
    }

    double sum_terms = 0.0;

    for (int j = left; j <= right; ++j) {
      double kernel_val =
        k_b_scalar(dist[i] - dist[j], b);

      sum_terms +=
        Z_v[j] *
        kernel_val *
        e[j] /
        (lambda * two_pi * dist[i]);
    }

    c0[i] = sum_terms / N_tau;
  }

  return c0;
}

