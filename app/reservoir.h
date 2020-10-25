#ifndef APP_RESERVOIR_H_
#define APP_RESERVOIR_H_

#define PRECISION_F32
//#define PRECISION_F64

#include "mat.h"

#if defined(PRECISION_F32)
#define VAL_T f32_t
#define MAT_T mat_f32_t
#define SPECTRAL_RADIUS_T float
#elif defined(PRECISION_F64)
#define VAL_T f64_t
#define MAT_T mat_f64_t
#define SPECTRAL_RADIUS_T double
#else
#error
#endif

#define DONT_TESET_XY 0
#define RESET_XY      1

typedef struct {
  mat_memory_t *mem;
  unsigned n_in_nodes;
  unsigned n_res_nodes;
  unsigned n_out_nodes;
  float leak_rate;
  MAT_T in_weights;  // heap: sizeof(VAL_T) * n_in_nodes * n_res_nodes
  MAT_T res_nodes;   // heap: sizeof(VAL_T) * 1 * n_res_nodes
  MAT_T res_weights; // heap: sizeof(VAL_T) * n_res_nodes * n_res_nodes
  MAT_T out_weights; // heap: sizeof(VAL_T) * n_res_nodes * n_out_nodes
  MAT_T x;           // heap: sizeof(VAL_T) * res->n_res_nodes * res->n_res_nodes
  MAT_T y;           // heap: sizeof(VAL_T) * res->n_in_nodes * res->n_res_nodes
} reservoir_t;

int init(reservoir_t *res);
void deinit(reservoir_t *res);

int train_feed_data(reservoir_t *res, MAT_T *data);
int train_compute_weight(reservoir_t *res, unsigned reset);
int train(reservoir_t *res, MAT_T *data, unsigned reset);

int predict(reservoir_t *res, MAT_T *predicted, MAT_T *data);

#endif /* APP_RESERVOIR_H_ */
