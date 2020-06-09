#ifndef RAINBOW_ASM_H
#define RAINBOW_ASM_H
#include <stdint.h>
#include "rainbow_blas.h"

void gf16mat_prod_512_32_ontheflybitsliced_asm(uint8_t *c, const uint8_t *matA, const uint8_t *b);
void gf16mat_prod_512_36_ontheflybitsliced_asm(uint8_t *c, const uint8_t *matA, const uint8_t *b);
void gf16mat_prod_16_32_ontheflybitsliced_asm(uint8_t *c, const uint8_t *matA, const uint8_t *b);
void gf16mat_prod_18_32_ontheflybitsliced_asm(uint8_t *c, const uint8_t *matA, const uint8_t *b);

#define gf16mat_prod_512_32_asm(c,a,b) gf16mat_prod_512_32_ontheflybitsliced_asm(c,a,b)
#define gf16mat_prod_512_36_asm(c,a,b) gf16mat_prod_512_36_ontheflybitsliced_asm(c,a,b)
#define gf16mat_prod_16_32_asm(c,a,b) gf16mat_prod_16_32_ontheflybitsliced_asm(c,a,b)
#define gf16mat_prod_18_32_asm(c,a,b) gf16mat_prod_18_32_ontheflybitsliced_asm(c,a,b)


//#define gf16mat_prod_512_32_asm(c,a,b) gf16mat_prod(c, a,512, 32, b)
//#define gf16mat_prod_512_36_asm(c,a,b) gf16mat_prod(c, a, 512, 36, b)
//#define gf16mat_prod_16_32_asm(c,a,b) gf16mat_prod(c, a, 16, 32, b)
//#define gf16mat_prod_18_32_asm(c,a,b) gf16mat_prod(c, a, 18, 32, b)

#endif  