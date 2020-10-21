#include "reservoir.h"

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

int init(reservoir_t *res) {
  if(MAT_NEW(res->mem, &res->in_weights, res->n_in_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &res->res_nodes, 1, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &res->res_weights, res->n_res_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &res->out_weights, res->n_res_nodes, res->n_out_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &res->x, res->n_res_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &res->y, res->n_res_nodes, res->n_in_nodes) < 0) {
    goto oom_fail;
  }
  // temporarily use res->x and res->y for work matrix
  MAT_T temp1;
  MAT_NEW(NULL, &temp1, 1, res->n_res_nodes);
  temp1.data = res->x.data;
  MAT_T temp2;
  MAT_NEW(NULL, &temp2, 1, res->n_res_nodes);
  temp2.data = res->y.data;
  // initialize in_weights
  for(unsigned n = 0; n < res->in_weights.n; n++) {
    for(unsigned m = 0; m < res->in_weights.m; m++) {
      *MAT(res->in_weights, n, m) = RANDOM_NORMAL(0.0, 1.0) < 0.0 ? -0.1 : 0.1;
    }
  }
  // initialize res_nodes
  MAT_ZEROS(&res->res_nodes);
  // initialize res_weights
  MAT_RANDOM_NORMAL(&res->res_weights, 0.0, 1.0);
  SPECTRAL_RADIUS_T spectral_radius = MAT_MAX_ABS_EIGENVAL(&res->res_weights, &temp1, &temp2, 100);
  MAT_MUL(&res->res_weights, &res->res_weights, 1.0 / spectral_radius);
  // out_weights
  MAT_ZEROS(&res->out_weights);
  // x
  MAT_ZEROS(&res->x);
  // y
  MAT_ZEROS(&res->y);
  return 0;
oom_fail:
  MAT_DESTROY(res->mem, &res->in_weights);
  MAT_DESTROY(res->mem, &res->res_nodes);
  MAT_DESTROY(res->mem, &res->res_weights);
  MAT_DESTROY(res->mem, &res->out_weights);
  MAT_DESTROY(res->mem, &res->x);
  MAT_DESTROY(res->mem, &res->y);
  return -1;
}

void deinit(reservoir_t *res) {
  MAT_DESTROY(res->mem, &res->in_weights);
  MAT_DESTROY(res->mem, &res->res_nodes);
  MAT_DESTROY(res->mem, &res->res_weights);
  MAT_DESTROY(res->mem, &res->out_weights);
  MAT_DESTROY(res->mem, &res->x);
  MAT_DESTROY(res->mem, &res->y);
}

static void get_next_node_state(reservoir_t *res, MAT_T *temp, MAT_T *next, MAT_T *curr, MAT_T *data) {
//  MAT_ZEROS(next);
//  MAT_PRODUCT(temp, data, &res->in_weights);
//  MAT_SUM(next, next, temp);
  MAT_PRODUCT(next, data, &res->in_weights);
  MAT_PRODUCT(temp, curr, &res->res_weights);
  MAT_SUM(next, next, temp);
  MAT_MUL(next, next, res->leak_rate);
  MAT_MUL(temp, curr, 1.0 - res->leak_rate);
  MAT_SUM(next, next, temp);
  for(unsigned n = 0; n < next->n; n++) {
    for(unsigned m = 0; m < next->m; m++) {
      *MAT(*next, n, m) = activate(*MAT(*next, n, m));
    }
  }
}

void train(reservoir_t *res, MAT_T *data) {
  MAT_T res_nodes, x, inv_x, y;
  if(MAT_NEW(res->mem, &res_nodes, data->n + 1, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &x, res->n_res_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &inv_x, res->n_res_nodes, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &y, res->n_res_nodes, res->n_in_nodes) < 0) {
    goto oom_fail;
  }
  // copy last res_nodes to initial res_nodes for training
  memcpy(res_nodes.data, res->res_nodes.data, sizeof(VAL_T) * res->n_res_nodes);
  // update res_nodes with given data for length of given data times
  // temporarily use x for work matrix
  MAT_T _temp;
  MAT_NEW(NULL, &_temp, 1, res->n_res_nodes);
  _temp.data = x.data;
  MAT_T _curr, _next, _data;
  MAT_NEW(NULL, &_curr, 1, res->n_res_nodes);
  MAT_NEW(NULL, &_next, 1, res->n_res_nodes);
  MAT_NEW(NULL, &_data, 1, res->n_in_nodes);
  for(unsigned n = 0; n < data->n; n++) {
    _curr.data = res_nodes.data + res->n_res_nodes * n;
    _next.data = res_nodes.data + res->n_res_nodes * (n + 1);
    _data.data = data->data + n * res->n_in_nodes;
    get_next_node_state(res, &_temp, &_next, &_curr, &_data);
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
  // compute out weights with accumulated (Y_TARGET X_T) and (X X_T)
  // temporarily use inv_x for (identity matrix * beta)
  MAT_NEW(NULL, &_temp, res->n_res_nodes, res->n_res_nodes);
  _temp.data = inv_x.data;
  MAT_IDENTITY(&_temp, 0.1);
  MAT_SUM(&x, &res->x, &_temp);
  MAT_INV(&inv_x, &x);
  MAT_PRODUCT(&res->out_weights, &inv_x, &res->y);
  // rewind one
  res_nodes.data -= res->n_res_nodes;
  res_nodes.n += 1;
  // copy back last res_nodes
  memcpy(res->res_nodes.data, res_nodes.data + res->n_res_nodes * data->n, sizeof(VAL_T) * res->n_res_nodes);
oom_fail:
  MAT_DESTROY(res->mem, &res_nodes);
  MAT_DESTROY(res->mem, &x);
  MAT_DESTROY(res->mem, &inv_x);
  MAT_DESTROY(res->mem, &y);
}

void predict(reservoir_t *res, MAT_T *predicted, MAT_T *data) {
  MAT_T temp, next;
  if(MAT_NEW(res->mem, &temp, 1, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  if(MAT_NEW(res->mem, &next, 1, res->n_res_nodes) < 0) {
    goto oom_fail;
  }
  get_next_node_state(res, &temp, &next, &res->res_nodes, data);
  MAT_COPY(&res->res_nodes, &next);
  MAT_PRODUCT(predicted, &next, &res->out_weights);
oom_fail:
  MAT_DESTROY(res->mem, &temp);
  MAT_DESTROY(res->mem, &next);
}
