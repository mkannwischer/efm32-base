#include "hal.h"
#include "em_device.h"
#include "em_cmu.h"
#include "em_gpio.h"
#include "em_usart.h"
#include "em_chip.h"


void hal_setup(const enum clock_mode clock)
{
    //TODO: clock setup
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

void hal_send_str(const char *in)
{
    int i;
    for (i = 0; in[i] != 0; i++)
    {
        USART_Tx(USART5, *(unsigned char *)(in + i));
    }
    USART_Tx(USART5, '\n');
}