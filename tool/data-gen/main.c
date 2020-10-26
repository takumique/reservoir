#include <stdio.h>
#include <stdlib.h>

void main() {
  FILE *fp;
  fp = fopen("training.txt", "r");
  if(fp == NULL) {
    return;
  }
  float *data = malloc(sizeof(float) * 960);
  if(data == NULL) {
    fclose(fp);
    return;
  }
  char buf[32];
  for(unsigned i = 0; fgets(buf, sizeof(buf), fp); i++) {
    data[i] = strtof(buf, NULL);
  }
  fclose(fp);
  fp = fopen("training_data.c", "w");
  if(fp == NULL) {
    free(data);
    return;
  }
  fprintf(fp, "#include \"reservoir.h\"\n");
  fprintf(fp, "const VAL_T __training_data[] = {\n");
  for(unsigned i = 0; i < 960; i++) {
    if(i > 0) {
      fprintf(fp, ",\n");
    }
    fprintf(fp, "  %f", data[i]);
  }
  fprintf(fp, "\n};\n");
  fclose(fp);
  free(data);
}
