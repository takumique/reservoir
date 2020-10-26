#include <stdio.h>
#include <stdlib.h>

#include "reservoir.h"

void main() {
  FILE *fp;
  mat_memory_t mem = {
      .memory_alloc = (void *(*)(unsigned)) malloc,
      .memory_free = (void (*)(void *)) free,
  };
  reservoir_t res = {
      .mem = &mem,
      .n_in_nodes = 1,
      .n_res_nodes = 100,
      .n_out_nodes = 1,
      .leak_rate = 0.02f,
  };
  init(&res);
  fp = fopen("weights.c", "w");
  fprintf(fp, "#include \"reservoir.h\"\n");
  // in_weights
  fprintf(fp, "const VAL_T __in_weights[] = {\n");
  for(unsigned i = 0; i < res.in_weights.n * res.in_weights.m; i++) {
    if(i > 0) {
      fprintf(fp, ",\n");
    }
    fprintf(fp, "  %f", *(res.in_weights.data + i));
  }
  fprintf(fp, "\n};\n");
  // res_weights
  fprintf(fp, "const VAL_T __res_weights[] = {\n");
  for(unsigned i = 0; i < res.res_weights.n * res.res_weights.m; i++) {
    if(i > 0) {
      fprintf(fp, ",\n");
    }
    fprintf(fp, "  %f", *(res.res_weights.data + i));
  }
  fprintf(fp, "\n};\n");
  fclose(fp);
}
