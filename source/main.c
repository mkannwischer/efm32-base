/**************************************************************************/ /**
 * @main_Series0.c
 * @brief This project reads input on the UART RX port and echoes it back to
 * the TX port one line at a time.  It does not use the VCOM or stdio.  See
 * readme.txt for details.
 * @version 0.0.1
 ******************************************************************************
 * @section License
 * <b>Copyright 2018 Silicon Labs, Inc. http://www.silabs.com</b>
 *******************************************************************************
 *
 * This file is licensed under the Silabs License Agreement. See the file
 * "Silabs_License_Agreement.txt" for details. Before using this software for
 * any purpose, you must agree to the terms of that agreement.
 *
 ******************************************************************************/
#include "em_device.h"
#include "em_cmu.h"
#include "em_gpio.h"
#include "em_usart.h"
#include "em_chip.h"`

#define BUFFER_SIZE 80
char buffer[BUFFER_SIZE];

static void hal_init(void)
{
  USART_InitAsync_TypeDef init = USART_INITASYNC_DEFAULT;
  // Chip errata
  CHIP_Init();
  // Enable oscillator to GPIO and USART1 modules
  CMU_ClockEnable(cmuClock_GPIO, true);
  CMU_ClockEnable(cmuClock_USART5, true);

  // set pin modes for UART TX and RX pins
  //GPIO_PinModeSet(gpioPortE, 9, gpioModeInput, 0);
  GPIO_PinModeSet(gpioPortE, 8, gpioModePushPull, 1);

  // Initialize USART asynchronous mode and route pins
  USART_InitAsync(USART5, &init);
  // USART5->ROUTEPEN |= USART_ROUTEPEN_TXPEN | USART_ROUTEPEN_RXPEN;
  USART5->ROUTEPEN |= USART_ROUTEPEN_TXPEN;
}

static void hal_send_str(unsigned char *in)
{
  int i;
  for (i = 0; in[i] != 0; i++)
  {
    USART_Tx(USART5, *(unsigned char *)(in + i));
  }
  USART_Tx(USART5, '\n');
}

int main(void)
{
  hal_init();
  hal_send_str("Hello. I'm a Giant Gecko and I would like to play!");

  while (1)
    ;
}