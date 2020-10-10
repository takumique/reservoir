#include "reservoir.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "linalg.h"

#ifdef SINGLE_PRECISION
float activate(float a) {
  return tanhf(a);
}
#else
double activate(double a) {
  return tanh(a);
}
#endif

static void init_in_weights(reservoir_t *res) {
  for(unsigned i = 0; i < res->n_in_nodes * res->n_res_nodes; i++) {
    *(res->in_weights + i) = RANDOM_NORMAL(0.0, 1.0) < 0.0 ? -0.1 : 0.1;
  }
}

static void init_res_nodes(reservoir_t *res) {
  memset(res->res_nodes, 0, sizeof(NODE_T) * res->n_res_nodes);
}

static void init_res_weights(reservoir_t *res) {
  WEIGHT_T *x, *y;
#ifndef STATIC_MEMORY
  x = (WEIGHT_T *) res->mem->memory_alloc(sizeof(WEIGHT_T) * res->n_res_nodes);
  if(!x) {
    goto oom_fail;
  }
  y = (WEIGHT_T *) res->mem->memory_alloc(sizeof(WEIGHT_T) * res->n_res_nodes);
  if(!y) {
    goto oom_fail;
  }
#else
  x = (WEIGHT_T *) res->mem->work;
  y = (WEIGHT_T *) STATIC_MEMORY_OFFSET(res->mem->work, WEIGHT_T, res->n_res_nodes);
#endif
  MATRIX_RANDOM_NORMAL(res->res_weights, res->n_res_nodes, res->n_res_nodes, 0.0, 1.0);
  WEIGHT_T spectral_radius = MAX_ABS_EIGENVAL(res->res_weights, x, y, res->n_res_nodes, 100);
  for(unsigned i = 0; i < res->n_res_nodes * res->n_res_nodes; i++) {
    *(res->res_weights + i) /= spectral_radius;
  }
#ifndef STATIC_MEMORY
oom_fail:
  res->mem->memory_free(x);
  res->mem->memory_free(y);
#endif
}

static void init_out_weights(reservoir_t *res) {
  memset(res->out_weights, 0, sizeof(WEIGHT_T) * res->n_res_nodes * res->n_out_nodes);
}

int init(reservoir_t *res) {
#ifndef STATIC_MEMORY
  res->in_weights = res->mem->memory_alloc(sizeof(WEIGHT_T) * res->n_in_nodes * res->n_res_nodes);
  if(!res->in_weights) {
    goto oom_fail;
  }
  res->res_nodes = res->mem->memory_alloc(sizeof(WEIGHT_T) * res->n_res_nodes);
  if(!res->res_nodes) {
    goto oom_fail;
  }
  res->res_weights = res->mem->memory_alloc(sizeof(WEIGHT_T) * res->n_res_nodes * res->n_res_nodes);
  if(!res->res_weights) {
    goto oom_fail;
  }
  res->out_weights = res->mem->memory_alloc(sizeof(WEIGHT_T) * res->n_res_nodes * res->n_out_nodes);
  if(!res->out_weights) {
    goto oom_fail;
  }
  res->x = res->mem->memory_alloc(sizeof(NODE_T) * res->n_res_nodes * res->n_res_nodes);
  if(!res->x) {
    goto oom_fail;
  }
  res->y = res->mem->memory_alloc(sizeof(NODE_T) * res->n_in_nodes * res->n_res_nodes);
  if(!res->y) {
    goto oom_fail;
  }
#endif
  init_in_weights(res);
  init_res_nodes(res);
  init_res_weights(res);
  init_out_weights(res);
  // test
  memset(res->x, 0, sizeof(NODE_T) * res->n_res_nodes * res->n_res_nodes);
  memset(res->y, 0, sizeof(NODE_T) * res->n_in_nodes * res->n_res_nodes);
  return 0;
#ifndef STATIC_MEMORY
oom_fail:
  res->mem->memory_free(res->in_weights);
  res->mem->memory_free(res->res_weights);
  res->mem->memory_free(res->out_weights);
  res->mem->memory_free(res->x);
  res->mem->memory_free(res->y);
  return -1;
#endif
}

void deinit(reservoir_t *res) {
#ifndef STATIC_MEMORY
  res->mem->memory_free(res->in_weights);
  res->mem->memory_free(res->res_weights);
  res->mem->memory_free(res->out_weights);
#endif
}

static void get_next_node_state(reservoir_t *res, NODE_T *next, NODE_T *curr, DATA_T *data) {
  NODE_T *temp;
#ifndef STATIC_MEMORY
  temp = (NODE_T *) res->mem->memory_alloc(sizeof(NODE_T) * res->n_res_nodes);
  if(!temp) {
    goto oom_fail;
  }
#else
  temp = (NODE_T *) res->mem->work;
#endif
  memset(next, 0, sizeof(NODE_T) * res->n_res_nodes);
  MATRIX_PRODUCT(temp, data, res->in_weights, 1, 1, res->n_res_nodes);
  for(unsigned i = 0; i < res->n_res_nodes; i++) {
    *(next + i) += *(temp + i);
  }
  MATRIX_PRODUCT(temp, curr, res->res_weights, 1, res->n_res_nodes, res->n_res_nodes);
  for(unsigned i = 0; i < res->n_res_nodes; i++) {
    *(next + i) += *(temp + i);
  }
  for(unsigned i = 0; i < res->n_res_nodes; i++) {
    *(next + i) *= (DATA_T) res->leak_rate;
  }
  VECTOR_MULTIPLICATION(temp, curr, 1.0 - res->leak_rate, res->n_res_nodes);
  for(unsigned i = 0; i < res->n_res_nodes; i++) {
    *(next + i) += *(temp + i);
  }
  for(unsigned i = 0; i < res->n_res_nodes; i++) {
    *(next + i) = activate(*(next + i));
  }
#ifndef STATIC_MEMORY
oom_fail:
  res->mem->memory_free(temp);
#endif
}

void train_batch(reservoir_t *res, DATA_T *data, unsigned data_len) {
  NODE_T *res_nodes_train, *x, *inv_x, *y, beta;
#ifndef STATIC_MEMORY
  res_nodes_train = (NODE_T *) res->mem->memory_alloc(sizeof(NODE_T) * res->n_res_nodes * (data_len + 1));
  if(!res_nodes_train) {
    goto oom_fail;
  }
  x = (NODE_T *) res->mem->memory_alloc(sizeof(NODE_T) * res->n_res_nodes * res->n_res_nodes);
  if(!x) {
    goto oom_fail;
  }
  inv_x = (NODE_T *) res->mem->memory_alloc(sizeof(WEIGHT_T) * res->n_res_nodes * res->n_res_nodes);
  if(!inv_x) {
    goto oom_fail;
  }
  y = (NODE_T *) res->mem->memory_alloc(sizeof(NODE_T) * res->n_in_nodes * res->n_res_nodes);
  if(!y) {
    goto oom_fail;
  }
#else
  res_nodes_train = (NODE_T *) res->mem->work;
  x = (NODE_T *) STATIC_MEMORY_OFFSET(res_nodes_train, NODE_T, res->n_res_nodes * (data_len + 1));
  inv_x = (NODE_T *) STATIC_MEMORY_OFFSET(x, NODE_T, res->n_res_nodes * n_res_nodes);
  y = (NODE_T *) STATIC_MEMORY_OFFSET(x, NODE_T, res->n_res_nodes * n_res_nodes);
#endif

#ifdef SINGLE_PRECISION
  beta = 0.1f;
#else
  beta = 0.1;
#endif

  // input loop
  memcpy(res_nodes_train, res->res_nodes, sizeof(NODE_T) * res->n_res_nodes);
  for(unsigned i = 0; i < data_len; i++) {
    NODE_T *prev = res_nodes_train + res->n_res_nodes * i;
    NODE_T *curr = res_nodes_train + res->n_res_nodes * (i + 1);
    get_next_node_state(res, curr, prev, data + i);
  }
  NODE_T *head = res_nodes_train + res->n_res_nodes;
  // update (X X_T)
  MATRIX_PRODUCT_LEFT_T(x, head, head, res->n_res_nodes, data_len, res->n_res_nodes);
  for(unsigned i = 0; i < res->n_res_nodes * res->n_res_nodes; i++) {
    *(res->x + i) += *(x + i);
  }
  memcpy(x, res->x, sizeof(NODE_T) * res->n_res_nodes * res->n_res_nodes);
  // update (Y_TARGET X_T)
  MATRIX_PRODUCT_LEFT_T(y, head, data, res->n_res_nodes, data_len, res->n_in_nodes);
  for(unsigned i = 0; i < res->n_in_nodes * res->n_res_nodes; i++) {
    *(res->y + i) += *(y + i);
  }
  // compute out weights with accumulated (Y_TARGET X_T) and (X X_T)
  for(unsigned i = 0; i < res->n_res_nodes; i++) {
    *(x + i * res->n_res_nodes + i) += beta;
  }
  INV(inv_x, x, res->n_res_nodes);
  MATRIX_PRODUCT(res->out_weights, inv_x, res->y, res->n_res_nodes, res->n_res_nodes, res->n_out_nodes);
  // copy reservoir nodes at the end of training
  memcpy(res->res_nodes, res_nodes_train + res->n_res_nodes * data_len, sizeof(NODE_T) * res->n_res_nodes);
#ifndef STATIC_MEMORY
oom_fail:
  res->mem->memory_free(res_nodes_train);
  res->mem->memory_free(x);
  res->mem->memory_free(inv_x);
  res->mem->memory_free(y);
#endif
}

void predict(reservoir_t *res, DATA_T *predicted, DATA_T data) {
  NODE_T *next;
#ifndef STATIC_MEMORY
  next = (NODE_T *) res->mem->memory_alloc(sizeof(NODE_T) * res->n_res_nodes);
  if(!next) {
    goto oom_fail;
  }
#else
  next = (NODE_T *) res->mem->work;
#endif
  get_next_node_state(res, next, res->res_nodes, &data);
  memcpy(res->res_nodes, next, sizeof(NODE_T) * res->n_res_nodes);
  MATRIX_PRODUCT(predicted, next, res->out_weights, 1, res->n_res_nodes, 1);
#ifndef STATIC_MEMORY
oom_fail:
  res->mem->memory_free(next);
#endif
}
