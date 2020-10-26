#include "reservoir.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG_PRINT

static inline void print_res_head(reservoir_t *res, char *label) {
#ifdef DEBUG_PRINT
  if(label) {
    printf("%s\n", label);
  }
  printf("res->in_weights: %f %f %f %f\n",
      *MAT(res->in_weights, 0, 0),
      *MAT(res->in_weights, 0, 1),
      *MAT(res->in_weights, 0, 2),
      *MAT(res->in_weights, 0, 3));
  printf("res->res_weights: %f %f %f %f\n",
      *MAT(res->res_weights, 0, 0),
      *MAT(res->res_weights, 0, 1),
      *MAT(res->res_weights, 0, 2),
      *MAT(res->res_weights, 0, 3));
  printf("res->out_weights: %f %f %f %f\n",
      *MAT(res->out_weights, 0, 0),
      *MAT(res->out_weights, 0, 1),
      *MAT(res->out_weights, 0, 2),
      *MAT(res->out_weights, 0, 3));
  printf("res->res_nodes: %f %f %f %f\n",
      *MAT(res->res_nodes, 0, 0),
      *MAT(res->res_nodes, 0, 1),
      *MAT(res->res_nodes, 0, 2),
      *MAT(res->res_nodes, 0, 3));
  printf("res->x: %f %f %f %f\n",
      *MAT(res->x, 0, 0),
      *MAT(res->x, 0, 1),
      *MAT(res->x, 0, 2),
      *MAT(res->x, 0, 3));
  printf("res->y: %f %f %f %f\n",
      *MAT(res->y, 0, 0),
      *MAT(res->y, 0, 1),
      *MAT(res->y, 0, 2),
      *MAT(res->y, 0, 3));
#endif
}

#define INPUT_FILE_NAME "training.txt"
#define INPUT_FILE_LINE_MAX 32
#define OUTPUT_FILE_NAME "predicted.txt"

#define TRAINING_DATA_SIZE 960
//#define TRAINING_BATCH_SIZE TRAINING_DATA_SIZE
#define TRAINING_BATCH_SIZE 10
#define PREDICTION_DATA_SIZE 640

int main(int argc, char** argv) {
  FILE *fp;
  mat_memory_t mem = {
      .memory_alloc = (void *(*)(unsigned)) malloc,
      .memory_free = (void (*)(void *)) free,
  };
  reservoir_t res = {
      .mem = &mem,
      .n_in_nodes = 1,
      .n_res_nodes = 150,
      .n_out_nodes = 1,
      .leak_rate = 0.02f,
  };
  // init
  init(&res);
  print_res_head(&res, "INIT");
  // load data + train
  MAT_T training_data;
  if(MAT_NEW(res.mem, &training_data, TRAINING_BATCH_SIZE, 1) < 0) {
    goto error;
  }
  char buf[INPUT_FILE_LINE_MAX];
  fp = fopen(INPUT_FILE_NAME, "r");
  if(fp == NULL) {
    goto error;
  }
  VAL_T data;
  for(unsigned i = 0; fgets(buf, sizeof(buf), fp);) {
    data = strtod(buf, NULL);
    *MAT(training_data, i, 0) = data;
    if(++i >= TRAINING_BATCH_SIZE) {
      train_feed_data(&res, &training_data);
      i = 0;
    }
  }
  train_compute_weight(&res, RESET_XY);
  fclose(fp);
  print_res_head(&res, "TRAINED");
  // predict + save
  MAT_T predicted_data;
  if(MAT_NEW(res.mem, &predicted_data, 1, 1) < 0) {
    goto error;
  }
  MAT_T prev_data;
  if(MAT_NEW(res.mem, &prev_data, 1, 1) < 0) {
    goto error;
  }
  *MAT(prev_data, 0, 0) = data;
  fp = fopen(OUTPUT_FILE_NAME, "w");
  if(fp == NULL) {
    goto error;
  }
  for(unsigned i = 0; i < PREDICTION_DATA_SIZE; i++) {
    predict(&res, &predicted_data, &prev_data);
    fprintf(fp, "%f\n", *MAT(predicted_data, 0, 0));
    *MAT(prev_data, 0, 0) = *MAT(predicted_data, 0, 0);
  }
  fclose(fp);
  // done
  deinit(&res);
  return 0;
error:
  MAT_DESTROY(res.mem, &training_data);
  MAT_DESTROY(res.mem, &predicted_data);
  MAT_DESTROY(res.mem, &prev_data);
  return -1;
}
