#ifndef APP_CMSIS_APPLICATION_H_
#define APP_CMSIS_APPLICATION_H_

#include "stm32l4xx_hal.h"

#include "application_log.h"

void application_init(UART_HandleTypeDef *uart);
void application_loop();
void application_tick();

#endif /* APP_CMSIS_APPLICATION_H_ */
