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
      (res->in_weights)[0],
      (res->in_weights)[1],
      (res->in_weights)[2],
      (res->in_weights)[3]);
  printf("res->res_weights: %f %f %f %f\n",
      (res->res_weights)[0],
      (res->res_weights)[1],
      (res->res_weights)[2],
      (res->res_weights)[3]);
  printf("res->out_weights: %f %f %f %f\n",
      (res->out_weights)[0],
      (res->out_weights)[1],
      (res->out_weights)[2],
      (res->out_weights)[3]);
  printf("res->res_nodes: %f %f %f %f\n",
      (res->res_nodes)[0],
      (res->res_nodes)[1],
      (res->res_nodes)[2],
      (res->res_nodes)[3]);
#endif
}

int main(int argc, char** argv) {
  FILE *fp;
  memory_ops_t mem = {
#ifndef STATIC_MEMORY
      .memory_alloc = (void *(*)(unsigned)) malloc,
      .memory_free = (void (*)(void *)) free,
#else
      // TODO
#endif
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
  // read file
  fp = fopen("training.txt", "r");
  if(fp == NULL) {
    return -1;
  }
  char buf[32];
  DATA_T training_data[960];
  unsigned i;
  for(i = 0; fgets(buf, sizeof(buf), fp); i++) {
#ifdef SINGLE_PRECISION
    training_data[i] = strtof(buf, NULL);
#else
    training_data[i] = strtod(buf, NULL);
#endif
  }
  i -= 1;
  fclose(fp);
  // train model
  train(&res, training_data, i + 1);
  print_res_head(&res, "TRAINED");
  // predict by pre-trained model and save to file
  fp = fopen("predicted.txt", "w");
  DATA_T data = training_data[i];
  for(i = 0; i < 640; i++) {
    predict(&res, &data, data);
    fprintf(fp, "%f\n", data);
  }
  fclose(fp);
  // done
  deinit(&res);
  return 0;
}
