#include "mat.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

//#define CHECK_ARGS

int mat_f32_new(mat_memory_t *mem, mat_f32_t *a, unsigned n, unsigned m) {
  a->t = 0;
  a->n = n;
  a->m = m;
  if(mem && mem->memory_alloc) {
    a->data = (f32_t *) mem->memory_alloc(sizeof(f32_t) * n * m);
    if(a->data == NULL) {
      return -1;
    }
  } else {
    a->data = NULL;
  }
  return 0;
}

void mat_f32_destroy(mat_memory_t *mem, mat_f32_t *a) {
  if(mem && mem->memory_free) {
    mem->memory_free(a->data);
    a->data = NULL;
  }
}

int mat_f32_copy(mat_f32_t *dst, mat_f32_t *src) {
  memcpy(dst->data, src->data, sizeof(f32_t) * dst->n * dst->m);
  return 0;
}

int mat_f32_zeros(mat_f32_t *a) {
  memset(a->data, 0, sizeof(f32_t) * a->n * a->m);
  return 0;
}

int mat_f32_transpose(mat_f32_t *at, mat_f32_t *a) {
#ifdef CHECK_ARGS
  if(at->n != a->m) {
    return -1;
  }
  if(at->m != a->n) {
    return -1;
  }
#endif
  arm_matrix_instance_f32 _at, _a;
  arm_mat_init_f32(&_at, at->n, at->m, at->data);
  arm_mat_init_f32(&_a, a->n, a->m, a->data);
  arm_mat_trans_f32(&_a, &_at);
  at->t = 0;
  at->n = a->m;
  at->m = a->n;
  return 0;
}

int mat_f32_sum(mat_f32_t *c, mat_f32_t *a, mat_f32_t *b) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != a->m) {
    return -1;
  }
  if(c->n != b->n) {
    return -1;
  }
  if(c->m != b->m) {
    return -1;
  }
#endif
  arm_matrix_instance_f32 _c, _a, _b;
  c->t = 0;
  arm_mat_init_f32(&_c, c->n, c->m, c->data);
  arm_mat_init_f32(&_a, a->n, a->m, a->data);
  arm_mat_init_f32(&_b, b->n, b->m, b->data);
  if(arm_mat_add_f32(&_a, &_b, &_c) < 0) {
    return -1;
  }
  return 0;
}

int mat_f32_add(mat_f32_t *c, mat_f32_t *a, float l) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != a->m) {
    return -1;
  }
#endif
  c->t = 0;
  for(unsigned i = 0; i < c->n * c->m; i++) {
    *(c->data + i) = *(a->data + i) + l;
  }
  return 0;
}

int mat_f32_product(mat_f32_t *c, mat_f32_t *a, mat_f32_t *b) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != b->m) {
    return -1;
  }
  if(a->m != b->n) {
    return -1;
  }
#endif
  arm_matrix_instance_f32 _c, _a, _b;
  c->t = 0;
  arm_mat_init_f32(&_c, c->n, c->m, c->data);
  arm_mat_init_f32(&_a, a->n, a->m, a->data);
  arm_mat_init_f32(&_b, b->n, b->m, b->data);
  if(arm_mat_mult_f32(&_a, &_b, &_c) < 0) {
    return -1;
  }
  return 0;
}

int mat_f32_mul(mat_f32_t *c, mat_f32_t *a, float l) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != a->m) {
    return -1;
  }
#endif
  arm_matrix_instance_f32 _c, _a;
  c->t = 0;
  arm_mat_init_f32(&_c, c->n, c->m, c->data);
  arm_mat_init_f32(&_a, a->n, a->m, a->data);
  if(arm_mat_scale_f32(&_a, l, &_c) < 0) {
    return -1;
  }
  return 0;
}

int mat_f32_identity(mat_f32_t *c, f32_t l) {
#ifdef CHECK_ARGS
  if(c->n != c->m) {
    return -1;
  }
#endif
  c->t = 0;
  mat_f32_zeros(c);
  for(unsigned i = 0; i < c->n; i++) {
    *_MAT(*c, i, i) = l;
  }
  return 0;
}

int mat_f32_inv(mat_f32_t *inv_a, mat_f32_t *a) {
#ifdef CHECK_ARGS
  if(inv_a->n != a->n) {
    return -1;
  }
  if(inv_a->m != a->m) {
    return -1;
  }
  if(inv_a->n != inv_a->m) {
    return -1;
  }
#endif
  arm_matrix_instance_f32 _inv_a, _a;
  inv_a->t = 0;
  arm_mat_init_f32(&_inv_a, inv_a->n, inv_a->m, inv_a->data);
  arm_mat_init_f32(&_a, a->n, a->m, a->data);
  if(arm_mat_inverse_f32(&_a, &_inv_a) < 0) {
    return -1;
  }
  return 0;
}

f32_t f32_random_normal(float mu, float sigma) {
  float z = sqrt(-2.0f * log((float) random() / RAND_MAX)) * cos(2.0f * M_PI * ((float) random() / RAND_MAX));
  return mu + sigma * z;
}

void mat_f32_random_normal(mat_f32_t *c, float mu, float sigma) {
  c->t = 0;
  for(unsigned n = 0; n < c->n; n++) {
    for(unsigned m = 0; m < c->m; m++) {
      *_MAT(*c, n, m) = f32_random_normal(mu, sigma);
    }
  }
}

float mat_f32_max_abs_eigenval(mat_f32_t *a, mat_f32_t *x, mat_f32_t *y, unsigned lim) {
#ifdef CHECK_ARGS
  if(a->n != a->m) {
    return -1;
  }
  if(x->n == 0 || a->m != x->m) {
    return -1;
  }
  if(y->n == 0 || a->m != y->m) {
    return -1;
  }
#endif
  float lambda = 1.0f, prev_lambda = 0.0f;
  for(unsigned i = 0; i < a->n; i++) {
    *_MAT(*x, 0, i) = 1.0f;
  }
  for(unsigned iteration = 0; iteration < lim; iteration++) {
    for(unsigned i = 0; i < a->n; i++) {
      *_MAT(*y, 0, i) = 0.0f;
      for(unsigned j = 0; j < a->n; j++) {
        *_MAT(*y, 0, i) += *_MAT(*a, i, j) * *_MAT(*x, 0, j);
      }
    }
    for(unsigned i = 0; i < a->n; i++){
      if(fabs(*_MAT(*y, 0, i)) > fabs(prev_lambda)) {
        lambda = *_MAT(*y, 0, i);
      }
    }
    for(unsigned i = 0; i < a->n; i++){
      *_MAT(*x, 0, i) = *_MAT(*y, 0, i) / lambda;
    }
    if(fabs(lambda - prev_lambda) < 10e-6f) {
      break;
    }
    prev_lambda = lambda;
  }
  return lambda;
}

int mat_f64_new(mat_memory_t *sup, mat_f64_t *a, unsigned n, unsigned m) {
  a->t = 0;
  a->n = n;
  a->m = m;
  if(sup && sup->memory_alloc) {
    a->data = (f64_t *) sup->memory_alloc(sizeof(f64_t) * n * m);
    if(a->data == NULL) {
      return -1;
    }
  } else {
    a->data = NULL;
  }
  return 0;
}

void mat_f64_destroy(mat_memory_t *mem, mat_f64_t *a) {
  if(mem && mem->memory_free) {
    mem->memory_free(a->data);
    a->data = NULL;
  }
}

int mat_f64_copy(mat_f64_t *dst, mat_f64_t *src) {
  memcpy(dst->data, src->data, sizeof(f32_t) * dst->n * dst->m);
  return 0;
}

int mat_f64_zeros(mat_f64_t *a) {
  memset(a->data, 0, sizeof(f64_t) * a->n * a->m);
  return 0;
}

int mat_f64_transpose(mat_f64_t *at, mat_f64_t *a) {
#ifdef CHECK_ARGS
  if(at->n != a->m) {
    return -1;
  }
  if(at->m != a->n) {
    return -1;
  }
#endif
  mat_f64_copy(at, a);
  at->t = !a->t;
  at->n = a->m;
  at->m = a->n;
  return 0;
}

int mat_f64_sum(mat_f64_t *c, mat_f64_t *a, mat_f64_t *b) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != a->m) {
    return -1;
  }
  if(c->n != b->n) {
    return -1;
  }
  if(c->m != b->m) {
    return -1;
  }
#endif
  for(;;) {
    if(c->n == 1 || c->m == 1) {
      for(unsigned i = 0; i < c->n * c->m; i++) {
        *(c->data + i) = *(a->data + i) + *(b->data + i);
      }
      break;
    }
    if(!a->t && !b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = *_MAT(*a, n, m) + *_MAT(*b, n, m);
        }
      }
      break;
    }
    if(a->t && !b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = *_MAT_T(*a, n, m) + *_MAT(*b, n, m);
        }
      }
      break;
    }
    if(!a->t && b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = *_MAT(*a, n, m) + *_MAT_T(*b, n, m);
        }
      }
      break;
    }
    if(a->t && b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = *_MAT_T(*a, n, m) + *_MAT_T(*b, n, m);
        }
      }
      break;
    }
  }
  return 0;
}

int mat_f64_add(mat_f64_t *c, mat_f64_t *a, double l) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != a->m) {
    return -1;
  }
#endif
  for(unsigned i = 0; i < c->n * c->m; i++) {
    *(c->data + i) = *(a->data + i) + l;
  }
  return 0;
}

int mat_f64_product(mat_f64_t *c, mat_f64_t *a, mat_f64_t *b) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != b->m) {
    return -1;
  }
  if(a->m != b->n) {
    return -1;
  }
#endif
  c->t = 0;
  for(;;) {
    if(!a->t && !b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = 0;
          for(unsigned p = 0; p < a->m; p++) {
            *_MAT(*c, n, m) += *_MAT(*a, n, p) * *_MAT(*b, p, m);
          }
        }
      }
      break;
    }
    if(a->t && !b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = 0;
          for(unsigned p = 0; p < a->m; p++) {
            *_MAT(*c, n, m) += *_MAT_T(*a, n, p) * *_MAT(*b, p, m);
          }
        }
      }
      break;
    }
    if(!a->t && b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = 0;
          for(unsigned p = 0; p < a->m; p++) {
            *_MAT(*c, n, m) += *_MAT(*a, n, p) * *_MAT_T(*b, p, m);
          }
        }
      }
      break;
    }
    if(a->t && b->t) {
      for(unsigned n = 0; n < c->n; n++) {
        for(unsigned m = 0; m < c->m; m++) {
          *_MAT(*c, n, m) = 0;
          for(unsigned p = 0; p < a->m; p++) {
            *_MAT(*c, n, m) += *_MAT_T(*a, n, p) * *_MAT_T(*b, p, m);
          }
        }
      }
      break;
    }
  }
  return 0;
}

int mat_f64_mul(mat_f64_t *c, mat_f64_t *a, double l) {
#ifdef CHECK_ARGS
  if(c->n != a->n) {
    return -1;
  }
  if(c->m != a->m) {
    return -1;
  }
#endif
  for(unsigned i = 0; i < c->n * c->m; i++) {
    *(c->data + i) = *(a->data + i) * l;
  }
  return 0;
}

int mat_f64_identity(mat_f64_t *c, double l) {
#ifdef CHECK_ARGS
  if(c->n != c->m) {
    return -1;
  }
#endif
  c->t = 0;
  mat_f64_zeros(c);
  for(unsigned i = 0; i < c->n; i++) {
    *_MAT(*c, i, i) = l;
  }
  return 0;
}

int mat_f64_inv(mat_f64_t *inv_a, mat_f64_t *a) {
#ifdef CHECK_ARGS
  if(inv_a->n != a->n) {
    return -1;
  }
  if(inv_a->m != a->m) {
    return -1;
  }
  if(inv_a->n != inv_a->m) {
    return -1;
  }
#endif
  double temp;
  mat_f64_identity(inv_a, 1.0f);
  if(a->t) {
    for(unsigned i = 0; i < inv_a->n; i++) {
      temp = 1.0 / *_MAT_T(*a, i, i);
      for(unsigned j = 0; j < inv_a->n; j++) {
        *_MAT_T(*a, i, j) *= temp;
        *_MAT(*inv_a, i, j) *= temp;
      }
      for(unsigned j = 0; j < inv_a->n; j++) {
        if(i != j) {
          temp = *_MAT_T(*a, j, i);
          for(unsigned k = 0; k < inv_a->n; k++) {
            *_MAT_T(*a, j, k) -= *_MAT_T(*a, i, k) * temp;
            *_MAT(*inv_a, j, k) -= *_MAT(*inv_a, i, k) * temp;
          }
        }
      }
    }
  } else {
    for(unsigned i = 0; i < inv_a->n; i++) {
      temp = 1.0 / *_MAT(*a, i, i);
      for(unsigned j = 0; j < inv_a->n; j++) {
        *_MAT(*a, i, j) *= temp;
        *_MAT(*inv_a, i, j) *= temp;
      }
      for(unsigned j = 0; j < inv_a->n; j++) {
        if(i != j) {
          temp = *_MAT(*a, j, i);
          for(unsigned k = 0; k < inv_a->n; k++) {
            *_MAT(*a, j, k) -= *_MAT(*a, i, k) * temp;
            *_MAT(*inv_a, j, k) -= *_MAT(*inv_a, i, k) * temp;
          }
        }
      }
    }
  }
  return 0;
}

f64_t f64_random_normal(double mu, double sigma) {
  double z = sqrt(-2.0f * log((double) random() / RAND_MAX)) * cos(2.0f * M_PI * ((double) random() / RAND_MAX));
  return mu + sigma * z;
}

void mat_f64_random_normal(mat_f64_t *c, double mu, double sigma) {
  c->t = 0;
  for(unsigned n = 0; n < c->n; n++) {
    for(unsigned m = 0; m < c->m; m++) {
      *_MAT(*c, n, m) = f64_random_normal(mu, sigma);
    }
  }
}

double mat_f64_max_abs_eigenval(mat_f64_t *a, mat_f64_t *x, mat_f64_t *y, unsigned lim) {
#ifdef CHECK_ARGS
  if(a->n != a->m) {
    return -1;
  }
  if(x->n == 0 || a->m != x->m) {
    return -1;
  }
  if(y->n == 0 || a->m != y->m) {
    return -1;
  }
#endif
  double lambda = 1.0f, prev_lambda = 0.0;
  for(unsigned i = 0; i < a->n; i++) {
    *_MAT(*x, 0, i) = 1.0;
  }
  for(unsigned iteration = 0; iteration < lim; iteration++) {
    for(unsigned i = 0; i < a->n; i++) {
      *_MAT(*y, 0, i) = 0.0f;
      for(unsigned j = 0; j < a->n; j++) {
        *_MAT(*y, 0, i) += *_MAT(*a, i, j) * *_MAT(*x, 0, j);
      }
    }
    for(unsigned i = 0; i < a->n; i++){
      if(fabs(*_MAT(*y, 0, i)) > fabs(prev_lambda)) {
        lambda = *_MAT(*y, 0, i);
      }
    }
    for(unsigned i = 0; i < a->n; i++){
      *_MAT(*x, 0, i) = *_MAT(*y, 0, i) / lambda;
    }
    if(fabs(lambda - prev_lambda) < 10e-6) {
      break;
    }
    prev_lambda = lambda;
  }
  return lambda;
}

#if 0
#include <stdio.h>
#include <stdlib.h>
void main() {
  mat_memory_t mat_support = {
      .memory_alloc = (void *(*)(unsigned)) malloc,
      .memory_free = (void (*)(void *)) free,
  };
  mat_f32_t a;
  if(mat_f32_new(&mat_support, &a, 2, 2) != 0) {
    printf("mat_f32_new failed\n");
    return;
  }
  *MAT(a, 0, 0) = 1;
  *MAT(a, 0, 1) = 2;
  *MAT(a, 1, 0) = 3;
  *MAT(a, 1, 1) = 4;
  mat_f32_t b;
  if(mat_f32_new(&mat_support, &b, 2, 3) != 0) {
    printf("mat_f32_new failed\n");
    return;
  }
  *MAT(b, 0, 0) = 5;
  *MAT(b, 0, 1) = 6;
  *MAT(b, 0, 2) = 7;
  *MAT(b, 1, 0) = 8;
  *MAT(b, 1, 1) = 9;
  *MAT(b, 1, 2) = 10;
  mat_f32_t c;
  if(mat_f32_new(&mat_support, &c, 2, 3) != 0) {
    printf("mat_f32_new failed\n");
    return;
  }
  mat_f32_t d;
  if(mat_f32_new(&mat_support, &d, 3, 2) != 0) {
    printf("mat_f32_new failed\n");
    return;
  }

  // c = a @ b
  printf("c = a @ b\n");
  if(mat_f32_product(&c, &a, &b) != 0) {
    printf("mat_f32_product failed\n");
    return;
  }
  // expected:
  // [[21 24 27]
  //  [47 54 61]]
  printf("%f %f %f \n%f %f %f \n",
      *MAT(c, 0, 0),
      *MAT(c, 0, 1),
      *MAT(c, 0, 2),
      *MAT(c, 1, 0),
      *MAT(c, 1, 1),
      *MAT(c, 1, 2));

  // c = a.T @ b
  printf("c = a.T @ b\n");
  T(a);
  if(mat_f32_product(&c, &a, &b) != 0) {
    printf("mat_f32_product failed\n");
    return;
  }
  T(a);
  // expected:
  // [[29 33 37]
  //  [42 48 54]]
  printf("%f %f %f \n%f %f %f \n",
      *MAT(c, 0, 0),
      *MAT(c, 0, 1),
      *MAT(c, 0, 2),
      *MAT(c, 1, 0),
      *MAT(c, 1, 1),
      *MAT(c, 1, 2));

  // d = b.T @ a
  printf("d = b.T @ a\n");
  T(b);
  if(mat_f32_product(&d, &b, &a) != 0) {
    printf("mat_f32_product failed\n");
    return;
  }
  // expected:
  // [[29 42]
  //  [33 48]
  //  [37 54]]
  printf("%f %f \n%f %f \n%f %f \n",
      *MAT(d, 0, 0),
      *MAT(d, 0, 1),
      *MAT(d, 1, 0),
      *MAT(d, 1, 1),
      *MAT(d, 2, 0),
      *MAT(d, 2, 1));
}
#endif
