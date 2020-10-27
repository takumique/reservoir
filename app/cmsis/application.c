#include "application.h"

#include <stdlib.h>

#include "reservoir.h"

extern const VAL_T __training_data[];

static void print_res_head(reservoir_t *res, char *label) {
  if(label) {
	  LOG("%s\r\n", label);
  }
  LOG("res->in_weights: %f %f %f %f\r\n",
      *MAT(res->in_weights, 0, 0),
      *MAT(res->in_weights, 0, 1),
      *MAT(res->in_weights, 0, 2),
      *MAT(res->in_weights, 0, 3));
  LOG("res->res_weights: %f %f %f %f\r\n",
      *MAT(res->res_weights, 0, 0),
      *MAT(res->res_weights, 0, 1),
      *MAT(res->res_weights, 0, 2),
      *MAT(res->res_weights, 0, 3));
  LOG("res->out_weights: %f %f %f %f\r\n",
      *MAT(res->out_weights, 0, 0),
      *MAT(res->out_weights, 0, 1),
      *MAT(res->out_weights, 0, 2),
      *MAT(res->out_weights, 0, 3));
  LOG("res->res_nodes: %f %f %f %f\r\n",
      *MAT(res->res_nodes, 0, 0),
      *MAT(res->res_nodes, 0, 1),
      *MAT(res->res_nodes, 0, 2),
      *MAT(res->res_nodes, 0, 3));
  LOG("res->x: %f %f %f %f\r\n",
      *MAT(res->x, 0, 0),
      *MAT(res->x, 0, 1),
      *MAT(res->x, 0, 2),
      *MAT(res->x, 0, 3));
  LOG("res->y: %f %f %f %f\r\n",
      *MAT(res->y, 0, 0),
      *MAT(res->y, 0, 1),
      *MAT(res->y, 0, 2),
      *MAT(res->y, 0, 3));
}

void application_init(UART_HandleTypeDef *uart) {
  LOG_INIT(uart);

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
  // init
  init(&res);
  print_res_head(&res, "INIT");
  // train
  MAT_T training_data;
  MAT_NEW(NULL, &training_data, 1, 1);
  for(unsigned i = 0; i < 960; i++) {
    training_data.data = __training_data + i;
    train_feed_data(&res, &training_data);
  }
  train_compute_weight(&res, RESET_XY);
  print_res_head(&res, "TRAINED");

  // for indicating activity (turn on green LED on board)
  HAL_GPIO_WritePin(GPIOA, GPIO_PIN_5, 1);
}

// application main loop
void application_loop() {
}

// application tick (10ms timer handler)
void application_tick() {
}
