#ifndef APP_CMSIS_APPLICATION_LOG_H_
#define APP_CMSIS_APPLICATION_LOG_H_

#include <string.h>
#include <stdio.h>

#include "stm32l4xx_hal.h"

#define LINE_MAX (64)

void _log_init(UART_HandleTypeDef *log_uart);
void _log(const char *buf);

#define LOG_INIT(UART) _log_init(UART)
#define LOG(...) \
do { \
  char _buf[LINE_MAX];\
  snprintf(_buf, LINE_MAX, __VA_ARGS__);\
  _log(_buf); \
} while(0)

#endif /* APP_CMSIS_APPLICATION_LOG_H_ */
