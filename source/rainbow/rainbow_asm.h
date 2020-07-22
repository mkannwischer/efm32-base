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

extern void batch_quad_trimat_eval_gf16_100_32_asm(unsigned char * y, const unsigned char * trimat, const unsigned char * x);
extern void batch_quad_trimat_eval_gf16_96_32_asm(unsigned char * y, const unsigned char * trimat, const unsigned char * x);
extern void batch_quad_trimat_eval_gf16_32_16_asm( unsigned char * y, const unsigned char * trimat, const unsigned char * x);
extern void batch_quad_trimat_eval_gf16_36_16_asm( unsigned char * y, const unsigned char * trimat, const unsigned char * x);


//#define gf16mat_prod_512_32_asm(c,a,b) gf16mat_prod(c, a,512, 32, b)
//#define gf16mat_prod_512_36_asm(c,a,b) gf16mat_prod(c, a, 512, 36, b)
//#define gf16mat_prod_16_32_asm(c,a,b) gf16mat_prod(c, a, 16, 32, b)
//#define gf16mat_prod_18_32_asm(c,a,b) gf16mat_prod(c, a, 18, 32, b)

#endif  