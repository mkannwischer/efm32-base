#include "hal.h"

int main(void)
{
  hal_setup(CLOCK_FAST);
  hal_send_str("Hello. I'm a Giant Gecko and I would like to play!");

  while (1)
    ;
}