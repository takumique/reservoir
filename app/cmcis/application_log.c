#include "application_log.h"

UART_HandleTypeDef *_log_uart = NULL;

void _log_init(UART_HandleTypeDef *log_uart) {
  _log_uart = log_uart;
}

void _log(const char *buf) {
  if(_log_uart != NULL) {
    HAL_UART_Transmit(_log_uart, (uint8_t *) buf, strlen(buf), 0xffff);
  }
}
