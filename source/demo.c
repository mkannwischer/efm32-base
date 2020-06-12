#include "api.h"
#include "randombytes.h"
#include "hal.h"

#include <string.h>

int main(void)
{
    unsigned char str[100];
    volatile uint64_t t0, t1;
    uint32_t num;
    hal_setup(CLOCK_FAST);

    hal_send_str("==========================");
    hal_send_str("Hello. I'm a Giant Gecko and I would like to play!\n\n");


    // cycle counting
    hal_send_str("I can count cycles up to 2^64, look:");
    t0 = hal_get_time();
    for(int i=0;i< (1<<16);i++){
      asm("nop");
    }
    t1 = hal_get_time();
    sprintf(str, "Count: %d\n", t1-t0);
    hal_send_str(str);

    // TRNG
    hal_send_str("I can also produce true random numbers:");
    randombytes(&num, 4);
    sprintf(str, "Today's lucky number is: %u\n", num);
    hal_send_str(str);



    hal_send_str("#");
    while(1);
    return 0;
}
