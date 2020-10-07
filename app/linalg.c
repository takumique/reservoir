#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "linalg.h"

#define OFFSET(A, N, M, n, m) (A + M * n + m)
#define OFFSET_T(A, N, M, n, m) (A + N * m + n)

void matrix_product_f32(float *c, float *a, float *b, unsigned n, unsigned m, unsigned p) {
  for(unsigned _n = 0; _n < n; _n++) {
    for(unsigned _p = 0; _p < p; _p++) {
      *(c + _n * p + _p) = 0;
      for(unsigned _m = 0; _m < m; _m++) {
        *OFFSET(c, n, p, _n, _p) += *OFFSET(a, n, m, _n, _m) * *OFFSET(b, m, p, _m, _p);
      }
    }
  }
}

void matrix_product_left_t_f32(float *c, float *a, float *b, unsigned n, unsigned m, unsigned p) {
  for(unsigned _n = 0; _n < n; _n++) {
    for(unsigned _p = 0; _p < p; _p++) {
      *(c + _n * p + _p) = 0;
      for(unsigned _m = 0; _m < m; _m++) {
        *OFFSET(c, n, p, _n, _p) += *OFFSET_T(a, n, m, _n, _m) * *OFFSET(b, m, p, _m, _p);
      }
    }
  }
}

void matrix_product_right_t_f32(float *c, float *a, float *b, unsigned n, unsigned m, unsigned p) {
  for(unsigned _n = 0; _n < n; _n++) {
    for(unsigned _p = 0; _p < p; _p++) {
      *(c + _n * p + _p) = 0;
      for(unsigned _m = 0; _m < m; _m++) {
        *OFFSET(c, n, p, _n, _p) += *OFFSET(a, n, m, _n, _m) * *OFFSET_T(b, m, p, _m, _p);
      }
    }
  }
}

void vector_multiplication_f32(float *c, float *a, float l, unsigned n) {
  for(unsigned i = 0; i < n; i++) {
    *(c + i) = *(a + i) * l;
  }
}

void matrix_multiplication_f32(float *c, float *a, float l, unsigned n, unsigned m) {
  vector_multiplication_f32(c, a, l, n * m);
}

void identity_multiplied_f32(float *c, float r, unsigned n) {
  memset(c, 0, sizeof(float) * n * n);
  for(unsigned _n = 0; _n < n; _n++) {
    *(c + _n * n + _n) = r;
  }
}

void identity_f32(float *c, unsigned n) {
  identity_multiplied_f32(c, 1.0f, n);
}

void inv_f32(float *inv_a, float *a, unsigned n) {
  double buf;
  identity_f32(inv_a, n);
  for(unsigned i = 0; i < n; i++) {
    buf = 1.0f / *OFFSET(a, n, n, i, i);
    for(unsigned j = 0; j < n; j++) {
      *OFFSET(a, n, n, i, j) *= buf;
      *OFFSET(inv_a, n, n, i, j) *= buf;
    }
    for(unsigned j = 0; j < n; j++) {
      if(i != j) {
        buf = *OFFSET(a, n, n, j, i);
        for(unsigned k = 0; k < n; k++) {
          *OFFSET(a, n, n, j, k) -= *OFFSET(a, n, n, i, k) * buf;
          *OFFSET(inv_a, n, n, j, k) -= *OFFSET(inv_a, n, n, i, k) * buf;
        }
      }
    }
  }
}

float random_normal_f32(float mu, float sigma) {
  float z = sqrt(-2.0f * log((float) random() / RAND_MAX)) * cos(2.0f * M_PI * ((float) random() / RAND_MAX));
  return mu + sigma * z;
}

void vector_random_normal_f32(float *a, unsigned n, float mu, float sigma) {
  for(unsigned i = 0; i < n; i++) {
    *(a + i) = random_normal_f32(mu, sigma);
  }
}

void matrix_random_normal_f32(float *a, unsigned n, unsigned m, float mu, float sigma) {
  vector_random_normal_f32(a, n * m, mu, sigma);
}

float max_abs_eigenval_f32(float *a, float *x, float *y, unsigned n, unsigned lim) {
  float lambda, prev_lambda = 0.0;
  for(unsigned i = 0; i < n; i++) {
    x[i] = 1.0;
  }
  for(unsigned iteration = 0; iteration < lim; iteration++) {
    for(unsigned i = 0; i < n; i++) {
      *(y + i) = 0.0;
      for(unsigned j = 0; j < n; j++) {
        *(y + i) += *(a + i * n + j) * *(x + j);
      }
    }
    for(unsigned i = 0; i < n; i++){
      if(fabs(*(y + i)) > fabs(prev_lambda)) {
        lambda = *(y + i);
      }
    }
    for(unsigned i = 0; i < n; i++){
      *(x + i) = *(y + i)/lambda;
    }
    if(fabs(lambda - prev_lambda) < 10e-6) {
      break;
    }
    prev_lambda = lambda;
  }
  return lambda;
}

void matrix_product_f64(double *c, double *a, double *b, unsigned n, unsigned m, unsigned p) {
  for(unsigned _n = 0; _n < n; _n++) {
    for(unsigned _p = 0; _p < p; _p++) {
      *(c + _n * p + _p) = 0;
      for(unsigned _m = 0; _m < m; _m++) {
        *OFFSET(c, n, p, _n, _p) += *OFFSET(a, n, m, _n, _m) * *OFFSET(b, m, p, _m, _p);
      }
    }
  }
}

void matrix_product_left_t_f64(double *c, double *a, double *b, unsigned n, unsigned m, unsigned p) {
  for(unsigned _n = 0; _n < n; _n++) {
    for(unsigned _p = 0; _p < p; _p++) {
      *(c + _n * p + _p) = 0;
      for(unsigned _m = 0; _m < m; _m++) {
        *OFFSET(c, n, p, _n, _p) += *OFFSET_T(a, n, m, _n, _m) * *OFFSET(b, m, p, _m, _p);
      }
    }
  }
}

void matrix_product_right_t_f64(double *c, double *a, double *b, unsigned n, unsigned m, unsigned p) {
  for(unsigned _n = 0; _n < n; _n++) {
    for(unsigned _p = 0; _p < p; _p++) {
      *(c + _n * p + _p) = 0;
      for(unsigned _m = 0; _m < m; _m++) {
        *OFFSET(c, n, p, _n, _p) += *OFFSET(a, n, m, _n, _m) * *OFFSET_T(b, m, p, _m, _p);
      }
    }
  }
}

void vector_multiplication_f64(double *c, double *a, double l, unsigned n) {
  for(unsigned i = 0; i < n; i++) {
    *(c + i) = *(a + i) * l;
  }
}

void matrix_multiplication_f64(double *c, double *a, double l, unsigned n, unsigned m) {
  vector_multiplication_f64(c, a, l, n * m);
}

void identity_multiplied_f64(double *c, double r, unsigned n) {
  memset(c, 0, sizeof(double) * n * n);
  for(unsigned _n = 0; _n < n; _n++) {
    *(c + _n * n + _n) = r;
  }
}

void identity_f64(double *c, unsigned n) {
  identity_multiplied_f64(c, 1.0, n);
}

void inv_f64(double *inv_a, double *a, unsigned n) {
  double buf;
  identity_f64(inv_a, n);
  for(unsigned i = 0; i < n; i++) {
    buf = 1.0 / *OFFSET(a, n, n, i, i);
    for(unsigned j = 0; j < n; j++) {
      *OFFSET(a, n, n, i, j) *= buf;
      *OFFSET(inv_a, n, n, i, j) *= buf;
    }
    for(unsigned j = 0; j < n; j++) {
      if(i != j) {
        buf = *OFFSET(a, n, n, j, i);
        for(unsigned k = 0; k < n; k++) {
          *OFFSET(a, n, n, j, k) -= *OFFSET(a, n, n, i, k) * buf;
          *OFFSET(inv_a, n, n, j, k) -= *OFFSET(inv_a, n, n, i, k) * buf;
        }
      }
    }
  }
}

double random_normal_f64(double mu, double sigma) {
  double z = sqrt(-2 * log((double) random() / RAND_MAX)) * cos(2 * M_PI * ((double) random() / RAND_MAX));
  return mu + sigma * z;
}

void vector_random_normal_f64(double *a, unsigned n, double mu, double sigma) {
  for(unsigned i = 0; i < n; i++) {
    *(a + i) = random_normal_f64(mu, sigma);
  }
}

void matrix_random_normal_f64(double *a, unsigned n, unsigned m, double mu, double sigma) {
  vector_random_normal_f64(a, n * m, mu, sigma);
}

double max_abs_eigenval_f64(double *a, double *x, double *y, unsigned n, unsigned lim) {
  double lambda, prev_lambda = 0.0;
  for(unsigned i = 0; i < n; i++) {
    x[i] = 1.0;
  }
  for(unsigned iteration = 0; iteration < lim; iteration++) {
    for(unsigned i = 0; i < n; i++) {
      *(y + i) = 0.0;
      for(unsigned j = 0; j < n; j++) {
        *(y + i) += *(a + i * n + j) * *(x + j);
      }
    }
    for(unsigned i = 0; i < n; i++){
      if(fabs(*(y + i)) > fabs(prev_lambda)) {
        lambda = *(y + i);
      }
    }
    for(unsigned i = 0; i < n; i++){
      *(x + i) = *(y + i)/lambda;
    }
    if(fabs(lambda - prev_lambda) < 10e-6) {
      break;
    }
    prev_lambda = lambda;
  }
  return lambda;
}

