#ifndef APP_RESERVOIR_H_
#define APP_RESERVOIR_H_

#include <inttypes.h>

#define SINGLE_PRECISION

#ifdef SINGLE_PRECISION
#define DATA_T float
#define NODE_T float
#define WEIGHT_T float
#else
#define DATA_T double
#define NODE_T double
#define WEIGHT_T double
#endif

//#define STATIC_MEMORY

#ifdef STATIC_MEMORY
#define STATIC_MEMORY_OFFSET(BASE, TYPE, N) ((char *) BASE + sizeof(TYPE) * N)
#endif

typedef struct {
#ifndef STATIC_MEMORY
  void *(*memory_alloc)(unsigned);
  void (*memory_free)(void *);
#else
  void *work;
#endif
} memory_ops_t;

typedef struct {
  memory_ops_t *mem;
  unsigned n_in_nodes;
  unsigned n_res_nodes;
  unsigned n_out_nodes;
  float leak_rate;
  WEIGHT_T *in_weights;  // sizeof(WEIGHT_T) * n_in_nodes * n_res_nodes
  NODE_T *res_nodes;     // sizeof(NODE_T) * n_res_nodes
  WEIGHT_T *res_weights; // sizeof(WEIGHT_T) * n_res_nodes * n_res_nodes
  WEIGHT_T *out_weights; // sizeof(WEIGHT_T) * n_res_nodes * n_out_nodes
  NODE_T *x;             // (X X_T): sizeof(NODE_T) * res->n_res_nodes * res->n_res_nodes
  NODE_T *y;           // (Y_TARGET X_T): sizeof(NODE_T) * res->n_in_nodes * res->n_res_nodes
} reservoir_t;

int init(reservoir_t *res);
void deinit(reservoir_t *res);

void train(reservoir_t *res, DATA_T *data, unsigned data_len);
void train_batch(reservoir_t *res, DATA_T *data, unsigned data_len);
void predict(reservoir_t *res, DATA_T *predicted, DATA_T data);

#endif /* APP_RESERVOIR_H_ */
