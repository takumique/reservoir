#include "reservoir.h"

#ifdef CONST_WEIGHTS
extern const VAL_T __in_weights[];
extern const VAL_T __res_weights[];
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#if defined(PRECISION_F32)
static inline float activate(float a) {
  return tanhf(a);
}
#elif defined(PRECISION_F64)
static inline double activate(double a) {
  return tanh(a);
}
#endif

static void _init_in_weights(reservoir_t *res) {
  for(unsigned n = 0; n < res->in_weights.n; n++) {
    for(unsigned m = 0; m < res->in_weights.m; m++) {
      *_MAT(res->in_weights, n, m) = RANDOM_NORMAL(0.0, 1.0) < 0.0 ? -0.1 : 0.1;
    }
  }
}

static void _init_res_weights(reservoir_t *res, MAT_T *temp1, MAT_T *temp2) {
  MAT_RANDOM_NORMAL(&res->res_weights, 0.0, 1.0);
  SPECTRAL_RADIUS_T spectral_radius = MAT_MAX_ABS_EIGENVAL(&res->res_weights, temp1, temp2, 100);
  MAT_MUL(&res->res_weights, &res->res_weights, 1.0 / spectral_radius);
}

static void _init_xy(reservoir_t *res) {
  MAT_ZEROS(&res->x);
  MAT_ZEROS(&res->y);
}

int init(reservoir_t *res) {
#ifndef CONST_WEIGHTS
  MAT_T temp1, temp2;
#endif
#ifndef CONST_WEIGHTS
  if(MAT_NEW(res->mem, &res->in_weights, res->n_in_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
#else
  MAT_NEW(NULL, &res->in_weights, res->n_in_nodes, res->n_res_nodes);
  res->in_weights.data = (VAL_T *) __in_weights;
#endif
  if(MAT_NEW(res->mem, &res->res_nodes, 1, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
#ifndef CONST_WEIGHTS
  if(MAT_NEW(res->mem, &res->res_weights, res->n_res_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
#else
  MAT_NEW(NULL, &res->res_weights, res->n_res_nodes, res->n_res_nodes);
  res->res_weights.data = (VAL_T *) __res_weights;
#endif
  if(MAT_NEW(res->mem, &res->out_weights, res->n_res_nodes, res->n_out_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &res->x, res->n_res_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &res->y, res->n_res_nodes, res->n_in_nodes) < 0) {
    goto oom_fail;
  }
  srandom(0);
  // initialize in_weights
#ifndef CONST_WEIGHTS
  _init_in_weights(res);
#endif
  // initialize res_nodes
  MAT_ZEROS(&res->res_nodes);
  // initialize res_weights
#ifndef CONST_WEIGHTS
  MAT_NEW(NULL, &temp1, 1, res->n_res_nodes);
  temp1.data = res->x.data;
  MAT_NEW(NULL, &temp2, 1, res->n_res_nodes);
  temp2.data = res->y.data;
  _init_res_weights(res, &temp1, &temp2);
#endif
  // out_weights
  MAT_ZEROS(&res->out_weights);
  // x and y
  _init_xy(res);
  return 0;
oom_fail:
#ifndef CONST_WEIGHTS
  MAT_DESTROY(res->mem, &res->in_weights);
#endif
  MAT_DESTROY(res->mem, &res->res_nodes);
#ifndef CONST_WEIGHTS
  MAT_DESTROY(res->mem, &res->res_weights);
#endif
  MAT_DESTROY(res->mem, &res->out_weights);
  MAT_DESTROY(res->mem, &res->x);
  MAT_DESTROY(res->mem, &res->y);
  return -1;
}

void deinit(reservoir_t *res) {
#ifndef CONST_WEIGHTS
  MAT_DESTROY(res->mem, &res->in_weights);
#endif
  MAT_DESTROY(res->mem, &res->res_nodes);
#ifndef CONST_WEIGHTS
  MAT_DESTROY(res->mem, &res->res_weights);
#endif
  MAT_DESTROY(res->mem, &res->out_weights);
  MAT_DESTROY(res->mem, &res->x);
  MAT_DESTROY(res->mem, &res->y);
}

static void _get_next_node_state(reservoir_t *res, MAT_T *temp, MAT_T *next, MAT_T *curr, MAT_T *data) {
  MAT_PRODUCT(next, data, &res->in_weights);
  MAT_PRODUCT(temp, curr, &res->res_weights);
  MAT_SUM(next, next, temp);
  MAT_MUL(next, next, res->leak_rate);
  MAT_MUL(temp, curr, 1.0 - res->leak_rate);
  MAT_SUM(next, next, temp);
  for(unsigned i = 0; i < next->m; i++) {
    *(next->data + i) = activate(*(next->data + i));
  }
}

// heap: n_res_nodes * n_data + n_res_nodes * n_res_nodes + n_res_nodes * n_in_nodes
int train_feed_data(reservoir_t *res, MAT_T *data) {
  int ret = 0;
  MAT_T res_nodes, x, y;
  if(MAT_NEW(res->mem, &res_nodes, data->n + 1, res->n_res_nodes) < 0) {
    ret = -1;
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &x, res->n_res_nodes, res->n_res_nodes) < 0) {
    ret = -1;
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &y, res->n_res_nodes, res->n_in_nodes) < 0) {
    ret = -1;
    goto oom_fail;
  }
  // copy last res_nodes to initial res_nodes for training
  memcpy(res_nodes.data, res->res_nodes.data, sizeof(VAL_T) * res->n_res_nodes);
  // update res_nodes with given data for length of given data times
  MAT_T _curr, _next, _data, _temp;
  MAT_NEW(NULL, &_curr, 1, res->n_res_nodes);
  MAT_NEW(NULL, &_next, 1, res->n_res_nodes);
  MAT_NEW(NULL, &_data, 1, res->n_in_nodes);
  MAT_NEW(NULL, &_temp, 1, res->n_res_nodes);
  _temp.data = x.data;
  for(unsigned n = 0; n < data->n; n++) {
    _curr.data = res_nodes.data + res->n_res_nodes * n;
    _next.data = res_nodes.data + res->n_res_nodes * (n + 1);
    _data.data = data->data + n * res->n_in_nodes;
    _get_next_node_state(res, &_temp, &_next, &_curr, &_data);
  }
  // fast forward one
  res_nodes.data += res->n_res_nodes;
  res_nodes.n -= 1;
  // update (X X_T) as x = (X_T X) in case of column major
  MAT_T _res_nodes_t;
  MAT_NEW(NULL, &_res_nodes_t, res_nodes.n, res_nodes.m);
  _res_nodes_t.data = res_nodes.data;
  T(_res_nodes_t);
  MAT_PRODUCT(&x, &_res_nodes_t, &res_nodes);
  MAT_SUM(&res->x, &res->x, &x);
  // update (Y_TARGET X_T) as y = (X_T Y_TARGET) in case of column major
  MAT_PRODUCT(&y, &_res_nodes_t, data);
  MAT_SUM(&res->y, &res->y, &y);
  // rewind one
  res_nodes.data -= res->n_res_nodes;
  res_nodes.n += 1;
  // copy back last res_nodes
  memcpy(res->res_nodes.data, res_nodes.data + res->n_res_nodes * data->n, sizeof(VAL_T) * res->n_res_nodes);
oom_fail:
  MAT_DESTROY(res->mem, &res_nodes);
  MAT_DESTROY(res->mem, &x);
  MAT_DESTROY(res->mem, &y);
  return ret;
}

// heap: n_res_nodes * n_res_nodes * 1 (reset x and y)
// heap: n_res_nodes * n_res_nodes * 2 (do not reset x and y)
int train_compute_weight(reservoir_t *res, unsigned reset) {
  int ret = 0;
  MAT_T x, inv_x;
  if(!reset) {
    if(MAT_NEW(res->mem, &x, res->n_res_nodes, res->n_res_nodes) < 0) {
      ret = -1;
      goto oom_fail;
    }
  }
  if(MAT_NEW(res->mem, &inv_x, res->n_res_nodes, res->n_res_nodes) < 0) {
    ret = -1;
    goto oom_fail;
  }
  // compute out_weights with accumulated (Y_TARGET X_T) and (X X_T)
  if(!reset) {
    for(unsigned n = 0; n < x.n; n++) {
      for(unsigned m = 0; m < x.m; m++) {
        *_MAT(x, n, m) = *_MAT(res->x, n, m);
        if(n == m) {
          *_MAT(x, n, m) += 0.1;
        }
      }
    }
    MAT_INV(&inv_x, &x);
  } else {
    for(unsigned i = 0; i < res->x.n; i++) {
      *_MAT(res->x, i, i) += 0.1;
    }
    MAT_INV(&inv_x, &res->x);
  }
  MAT_PRODUCT(&res->out_weights, &inv_x, &res->y);
oom_fail:
  if(!reset) {
    MAT_DESTROY(res->mem, &x);
  } else {
    _init_xy(res);
  }
  MAT_DESTROY(res->mem, &inv_x);
  return ret;
}

int train(reservoir_t *res, MAT_T *data, unsigned reset) {
  if(train_feed_data(res, data) < 0) {
    return -1;
  }
  if(train_compute_weight(res, reset) < 0) {
    return -1;
  }
  return 0;
}

// heap: n_res_nodes + n_res_nodes
int predict(reservoir_t *res, MAT_T *predicted, MAT_T *data) {
  int ret = 0;
  MAT_T next, temp;
  if(MAT_NEW(res->mem, &next, 1, res->n_res_nodes) < 0) {
    ret = -1;
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &temp, 1, res->n_res_nodes) < 0) {
    ret = -1;
    goto oom_fail;
  }
  _get_next_node_state(res, &temp, &next, &res->res_nodes, data);
  MAT_COPY(&res->res_nodes, &next);
  MAT_PRODUCT(predicted, &res->res_nodes, &res->out_weights);
oom_fail:
  MAT_DESTROY(res->mem, &next);
  MAT_DESTROY(res->mem, &temp);
  return ret;
}
