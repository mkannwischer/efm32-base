#ifndef HAL_CONFIG_BOARD_H
#define HAL_CONFIG_BOARD_H

#include "em_device.h"
#include "hal-config-types.h"

// This file is auto-generated by Hardware Configurator in Simplicity Studio.
// Any content between $[ and ]$ will be replaced whenever the file is regenerated.
// Content outside these regions will be preserved.

// $[ACMP0]
// [ACMP0]$

// $[ACMP1]
// [ACMP1]$

// $[ANTDIV]
// [ANTDIV]$

// $[BTL_BUTTON]

#define BSP_BTL_BUTTON_PIN                            (2U)
#define BSP_BTL_BUTTON_PORT                           (gpioPortC)

// [BTL_BUTTON]$

// $[BUTTON]
#define BSP_BUTTON_PRESENT                            (1)

#define BSP_BUTTON0_PIN                               (2U)
#define BSP_BUTTON0_PORT                              (gpioPortC)

#define BSP_BUTTON1_PIN                               (3U)
#define BSP_BUTTON1_PORT                              (gpioPortC)

#define BSP_BUTTON_COUNT                              (2U)
#define BSP_BUTTON_INIT                               { { BSP_BUTTON0_PORT, BSP_BUTTON0_PIN }, { BSP_BUTTON1_PORT, BSP_BUTTON1_PIN } }
#define BSP_BUTTON_GPIO_DOUT                          (HAL_GPIO_DOUT_LOW)
#define BSP_BUTTON_GPIO_MODE                          (HAL_GPIO_MODE_INPUT)
// [BUTTON]$

// $[CMU]
#define BSP_CLK_HFXO_PRESENT                          (1)
#define BSP_CLK_HFXO_FREQ                             (38400000UL)
#define BSP_CLK_HFXO_INIT                              CMU_HFXOINIT_DEFAULT
#define BSP_CLK_HFXO_CTUNE                            (133)
#define BSP_CLK_LFXO_PRESENT                          (0)
#define BSP_CLK_LFXO_INIT                              CMU_LFXOINIT_DEFAULT
#define BSP_CLK_LFXO_FREQ                             (32768U)
#define BSP_CLK_LFXO_CTUNE                            (0U)
// [CMU]$

// $[COEX]
// [COEX]$

// $[EMU]
// [EMU]$

// $[EXTFLASH]
#define BSP_EXTFLASH_CS_PIN                           (3U)
#define BSP_EXTFLASH_CS_PORT                          (gpioPortD)

#define BSP_EXTFLASH_INTERNAL                         (0)
#define BSP_EXTFLASH_USART                            (HAL_SPI_PORT_USART0)
#define BSP_EXTFLASH_MOSI_PIN                         (0U)
#define BSP_EXTFLASH_MOSI_PORT                        (gpioPortD)

#define BSP_EXTFLASH_MISO_PIN                         (1U)
#define BSP_EXTFLASH_MISO_PORT                        (gpioPortD)

#define BSP_EXTFLASH_CLK_PIN                          (2U)
#define BSP_EXTFLASH_CLK_PORT                         (gpioPortD)

// [EXTFLASH]$

// $[EZRADIOPRO]
// [EZRADIOPRO]$

// $[GPIO]
#define PORTIO_GPIO_SWV_PIN                           (3U)
#define PORTIO_GPIO_SWV_PORT                          (gpioPortA)

#define BSP_TRACE_SWO_PIN                             (3U)
#define BSP_TRACE_SWO_PORT                            (gpioPortA)

// [GPIO]$

// $[I2C0]
#define PORTIO_I2C0_SCL_PIN                           (4U)
#define PORTIO_I2C0_SCL_PORT                          (gpioPortC)

#define PORTIO_I2C0_SDA_PIN                           (5U)
#define PORTIO_I2C0_SDA_PORT                          (gpioPortC)

#define BSP_I2C0_SCL_PIN                              (4U)
#define BSP_I2C0_SCL_PORT                             (gpioPortC)

#define BSP_I2C0_SDA_PIN                              (5U)
#define BSP_I2C0_SDA_PORT                             (gpioPortC)

// [I2C0]$

// $[I2C1]
// [I2C1]$

// $[I2CSENSOR]

#define BSP_I2CSENSOR_ENABLE_PIN                      (1U)
#define BSP_I2CSENSOR_ENABLE_PORT                     (gpioPortB)

#define BSP_I2CSENSOR_PERIPHERAL                      (HAL_I2C_PORT_I2C0)
#define BSP_I2CSENSOR_SCL_PIN                         (4U)
#define BSP_I2CSENSOR_SCL_PORT                        (gpioPortC)

#define BSP_I2CSENSOR_SDA_PIN                         (5U)
#define BSP_I2CSENSOR_SDA_PORT                        (gpioPortC)

// [I2CSENSOR]$

// $[IADC0]
// [IADC0]$

// $[IOEXP]
// [IOEXP]$

// $[LED]
#define BSP_LED_PRESENT                               (1)

#define BSP_LED0_PIN                                  (0U)
#define BSP_LED0_PORT                                 (gpioPortC)

#define BSP_LED1_PIN                                  (1U)
#define BSP_LED1_PORT                                 (gpioPortC)

#define BSP_LED_COUNT                                 (2U)
#define BSP_LED_INIT                                  { { BSP_LED0_PORT, BSP_LED0_PIN }, { BSP_LED1_PORT, BSP_LED1_PIN } }
#define BSP_LED_POLARITY                              (1)
// [LED]$

// $[LETIMER0]
// [LETIMER0]$

// $[LFXO]
// [LFXO]$

// $[MODEM]
// [MODEM]$

// $[PA]

#define BSP_PA_VOLTAGE                                (3300U)
// [PA]$

// $[PORTIO]
// [PORTIO]$

// $[PRS]
#define PORTIO_PRS_ASYNCH4_PIN                        (0U)
#define PORTIO_PRS_ASYNCH4_PORT                       (gpioPortB)

// [PRS]$

// $[PTI]
#define PORTIO_PTI_DFRAME_PIN                         (5U)
#define PORTIO_PTI_DFRAME_PORT                        (gpioPortC)

#define PORTIO_PTI_DOUT_PIN                           (4U)
#define PORTIO_PTI_DOUT_PORT                          (gpioPortC)

#define BSP_PTI_DFRAME_PIN                            (5U)
#define BSP_PTI_DFRAME_PORT                           (gpioPortC)

#define BSP_PTI_DOUT_PIN                              (4U)
#define BSP_PTI_DOUT_PORT                             (gpioPortC)

// [PTI]$

// $[SERIAL]
#define BSP_SERIAL_APP_PORT                           (HAL_SERIAL_PORT_USART1)
#define BSP_SERIAL_APP_TX_PIN                         (5U)
#define BSP_SERIAL_APP_TX_PORT                        (gpioPortA)

#define BSP_SERIAL_APP_RX_PIN                         (6U)
#define BSP_SERIAL_APP_RX_PORT                        (gpioPortA)

// [SERIAL]$

// $[SPIDISPLAY]

#define BSP_SPIDISPLAY_CS_PIN                         (4U)
#define BSP_SPIDISPLAY_CS_PORT                        (gpioPortD)

#define BSP_SPIDISPLAY_ENABLE_PIN                     (1U)
#define BSP_SPIDISPLAY_ENABLE_PORT                    (gpioPortB)

#define BSP_SPIDISPLAY_EXTCOMIN_PIN                   (0U)
#define BSP_SPIDISPLAY_EXTCOMIN_PORT                  (gpioPortB)

#define BSP_SPIDISPLAY_DISPLAY                        (HAL_DISPLAY_SHARP_LS013B7DH03)
#define BSP_SPIDISPLAY_USART                          (HAL_SPI_PORT_USART2)
#define BSP_SPIDISPLAY_EXTCOMIN_CHANNEL               (4)
#define BSP_SPIDISPLAY_MOSI_PIN                       (0U)
#define BSP_SPIDISPLAY_MOSI_PORT                      (gpioPortD)

#define BSP_SPIDISPLAY_MISO_PIN                       (1U)
#define BSP_SPIDISPLAY_MISO_PORT                      (gpioPortD)

#define BSP_SPIDISPLAY_CLK_PIN                        (2U)
#define BSP_SPIDISPLAY_CLK_PORT                       (gpioPortD)

// [SPIDISPLAY]$

// $[SPINCP]
#define BSP_SPINCP_NHOSTINT_PIN                       (2U)
#define BSP_SPINCP_NHOSTINT_PORT                      (gpioPortC)

#define BSP_SPINCP_NWAKE_PIN                          (3U)
#define BSP_SPINCP_NWAKE_PORT                         (gpioPortC)

#define BSP_SPINCP_USART_PORT                         (HAL_SPI_PORT_USART2)
#define BSP_SPINCP_MOSI_PIN                           (0U)
#define BSP_SPINCP_MOSI_PORT                          (gpioPortD)

#define BSP_SPINCP_MISO_PIN                           (1U)
#define BSP_SPINCP_MISO_PORT                          (gpioPortD)

#define BSP_SPINCP_CLK_PIN                            (2U)
#define BSP_SPINCP_CLK_PORT                           (gpioPortD)

#define BSP_SPINCP_CS_PIN                             (4U)
#define BSP_SPINCP_CS_PORT                            (gpioPortD)

// [SPINCP]$

// $[TIMER0]
// [TIMER0]$

// $[TIMER1]
// [TIMER1]$

// $[TIMER2]
// [TIMER2]$

// $[TIMER3]
// [TIMER3]$

// $[UARTNCP]
#define BSP_UARTNCP_USART_PORT                        (HAL_SERIAL_PORT_USART1)
#define BSP_UARTNCP_TX_PIN                            (5U)
#define BSP_UARTNCP_TX_PORT                           (gpioPortA)

#define BSP_UARTNCP_RX_PIN                            (6U)
#define BSP_UARTNCP_RX_PORT                           (gpioPortA)

// [UARTNCP]$

// $[USART0]
#define PORTIO_USART0_CLK_PIN                         (2U)
#define PORTIO_USART0_CLK_PORT                        (gpioPortD)

#define PORTIO_USART0_CS_PIN                          (3U)
#define PORTIO_USART0_CS_PORT                         (gpioPortD)

#define PORTIO_USART0_RX_PIN                          (1U)
#define PORTIO_USART0_RX_PORT                         (gpioPortD)

#define PORTIO_USART0_TX_PIN                          (0U)
#define PORTIO_USART0_TX_PORT                         (gpioPortD)

#define BSP_USART0_MOSI_PIN                           (0U)
#define BSP_USART0_MOSI_PORT                          (gpioPortD)

#define BSP_USART0_MISO_PIN                           (1U)
#define BSP_USART0_MISO_PORT                          (gpioPortD)

#define BSP_USART0_CLK_PIN                            (2U)
#define BSP_USART0_CLK_PORT                           (gpioPortD)

#define BSP_USART0_CS_PIN                             (3U)
#define BSP_USART0_CS_PORT                            (gpioPortD)

// [USART0]$

// $[USART1]
#define PORTIO_USART1_RX_PIN                          (6U)
#define PORTIO_USART1_RX_PORT                         (gpioPortA)

#define PORTIO_USART1_TX_PIN                          (5U)
#define PORTIO_USART1_TX_PORT                         (gpioPortA)

#define BSP_USART1_TX_PIN                             (5U)
#define BSP_USART1_TX_PORT                            (gpioPortA)

#define BSP_USART1_RX_PIN                             (6U)
#define BSP_USART1_RX_PORT                            (gpioPortA)

// [USART1]$

// $[USART2]
#define PORTIO_USART2_CLK_PIN                         (2U)
#define PORTIO_USART2_CLK_PORT                        (gpioPortD)

#define PORTIO_USART2_CS_PIN                          (4U)
#define PORTIO_USART2_CS_PORT                         (gpioPortD)

#define PORTIO_USART2_RX_PIN                          (1U)
#define PORTIO_USART2_RX_PORT                         (gpioPortD)

#define PORTIO_USART2_TX_PIN                          (0U)
#define PORTIO_USART2_TX_PORT                         (gpioPortD)

#define BSP_USART2_MOSI_PIN                           (0U)
#define BSP_USART2_MOSI_PORT                          (gpioPortD)

#define BSP_USART2_MISO_PIN                           (1U)
#define BSP_USART2_MISO_PORT                          (gpioPortD)

#define BSP_USART2_CLK_PIN                            (2U)
#define BSP_USART2_CLK_PORT                           (gpioPortD)

#define BSP_USART2_CS_PIN                             (4U)
#define BSP_USART2_CS_PORT                            (gpioPortD)

// [USART2]$

// $[VCOM]
// [VCOM]$

// $[VUART]
// [VUART]$

// $[WDOG]
// [WDOG]$

#if defined(_SILICON_LABS_MODULE)
#include "sl_module.h"
#endif

#endif /* HAL_CONFIG_BOARD_H */
