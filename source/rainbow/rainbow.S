.syntax unified
.cpu cortex-m4
.thumb


.macro mul2 a, tmp
  and \tmp, \a, #0x55555555
  and \a, \a, #0xaaaaaaaa
  eor \a, \a, \a, LSR#1
  eor \a, \a, \tmp, LSL#1
.endm

.global gf4v_mul_2_u32_asm
.type gf4v_mul_2_u32_asm, %function
.align 2
gf4v_mul_2_u32_asm:
  //uint32_t bit0 = a & 0x55555555;
  //uint32_t bit1 = a & 0xaaaaaaaa;
  //return (bit0 << 1) ^ bit1 ^ (bit1 >> 1);
  mul2 r0, r1
  bx lr

.global gf4v_mul_u32_asm
.type gf4v_mul_u32_asm, %function
.align 2
gf4v_mul_u32_asm:
    //uint32_t bit0_b = ((uint32_t) 0) - ((uint32_t)(b & 1));
    //uint32_t bit1_b = ((uint32_t) 0) - ((uint32_t)((b >> 1) & 1));
    //return (a & bit0_b) ^ (bit1_b & gf4v_mul_2_u32(a));

    sbfx r2, r1, #1, #1 // bit1_b
    sbfx r1, r1, #0, #1 // bit0_b
    and r1, r0, r1  // (a & bit0_b)

    //gf4v_mul_2_u32_asm
    mul2 r0, r3
    and r0, r2, r0 // bit1_b & gf4v_mul_2_u32(a)
    eor r0, r0, r1
    bx lr


.macro mulu32 a, b0, b1, b2, b3, t0, t1, t2
    and \t0, \a, \b0  // a&b0
    and \t1, \a, \b2  // a&b2
    mul2 \a, \t2

    and \t2, \a, \b1 // 2a & b1
    and \a, \a, \b3  // 2a & b3
    eor \t0, \t0, \t2
    eor \a, \t1, \a

    lsl \t1, \a, #2  // axb1
    and \t1, \t1, #0xcccccccc //a0b1
    eor \t0, \t0, \t1
    and \a, \a, #0xcccccccc //a1b1
    eor \t0, \t0, \a

    lsr \a, \a, #2 // a1b1_2

    //gf4v_mul_2_u32
    mul2 \a, \t1
    eor \a, \a, \t0
.endm

.macro accmulu32 a, c, b0, b1, b2, b3, t1, t2
    and \t2, \a, \b0  // a&b0
    eor \c, \c, \t2
    and \t1, \a, \b2  // a&b2
    mul2 \a, \t2

    and \t2, \a, \b1 // 2a & b1
    and \a, \a, \b3  // 2a & b3
    eor \c, \c, \t2
    eor \a, \t1, \a

    lsl \t1, \a, #2  // axb1
    and \t1, \t1, #0xcccccccc //a0b1
    eor \c, \c, \t1
    and \a, \a, #0xcccccccc //a1b1
    eor \c, \c, \a

    lsr \a, \a, #2 // a1b1_2

    //gf4v_mul_2_u32
    mul2 \a, \t1
    eor \c, \c, \a
.endm

.global gf16v_mul_u32_asm
.type gf16v_mul_u32_asm, %function
.align 2
gf16v_mul_u32_asm:
//uint32_t gf16v_mul_u32(uint32_t a, uint8_t b) {
//   uint32_t axb0 = gf4v_mul_u32(a, b);
//    uint32_t axb1 = gf4v_mul_u32(a, b >> 2);
//    uint32_t a0b1 = (axb1 << 2) & 0xcccccccc;
//    uint32_t a1b1 = axb1 & 0xcccccccc;
//    uint32_t a1b1_2 = a1b1 >> 2;
//    return axb0 ^ a0b1 ^ a1b1 ^ gf4v_mul_2_u32(a1b1_2);
//}
    push {r4-r12, lr}

    sbfx r4, r1, #3, #1  //b3
    sbfx r3, r1, #2, #1  //b2
    sbfx r2, r1, #1, #1  //b1
    sbfx r1, r1, #0, #1  //b0

    mulu32 r0, r1, r2, r3, r4, r10, r11, r12
    pop {r4-r12, lr}
    bx lr





.global gf16v_madd_u32_512_asm
.type gf16v_madd_u32_512_asm, %function
.align 2
gf16v_madd_u32_512_asm:
    push {r4-r11, r14}
    // r0 = accu
    // r1 = a
    // r2 = b

    sbfx r11, r2, #3, #1  //b3
    sbfx r10, r2, #2, #1  //b2
    sbfx r9, r2, #1, #1  //b1
    sbfx r8, r2, #0, #1  //b0

    // TODO: these values are wrong! it's 512 bytes
    // TODO: implement this as a loop
    // 42x3
    .rept 42
    ldmia r1!, {r2-r4}
    ldm r0, {r5-r7}
    accmulu32 r2, r5, r8, r9, r10, r11, r12, r14
    accmulu32 r3, r6, r8, r9, r10, r11, r12, r14
    accmulu32 r4, r7, r8, r9, r10, r11, r12, r14
    stmia r0!, {r5-r7}
    .endr
    // 1x2
    ldmia r1!, {r2-r3}
    ldm r0, {r5-r6}
    accmulu32 r2, r5, r8, r9, r10, r11, r12, r14
    accmulu32 r3, r6, r8, r9, r10, r11, r12, r14
    stmia r0!, {r5-r6}

    pop {r4-r11, r14}
    bx lr

//void gf16mat_prod_512_32_asm(uint8_t *c, const uint8_t *matA, const uint8_t *b);
.global gf16mat_prod_512_32_asm
.type gf16mat_prod_512_32_asm, %function
.align 2
gf16mat_prod_512_32_asm:

    // r0 = accu
    // r1 = A
    // r2 = b
    push {r4-r11, r14}



    mov r3, #0
    .set k, 0
    .rept 128
    str r3, [r0, k]
    .set k, k+4
    .endr

    one         .req s2
    ctr1        .req s3
    ctr2        .req s4

    vmov one, #1
    vmov ctr1, #4
    1:
    ldr r3, [r2], #4      // b0 ... b31
    vmov ctr2, #8
    2:

    sbfx r11, r3, #3, #1  // b3
    sbfx r10, r3, #2, #1  // b2
    sbfx r9, r3,  #1, #1  // b1
    sbfx r8, r3,  #0, #1  // b0
    .rept 64
    ldmia r1!, {r4-r5}
    ldm r0, {r6-r7}
    accmulu32 r4, r6, r8, r9, r10, r11, r12, r14
    accmulu32 r5, r7, r8, r9, r10, r11, r12, r14
    stmia r0!, {r6-r7}
    .endr

    sub r0, r0, #512
    lsr r3, r3,  #4
    vsub.f32 ctr2, ctr2, one
    vcmp.f32 ctr2, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 2b

    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    pop {r4-r11, r14}
    bx lr

    .unreq one
    .unreq ctr1
    .unreq ctr2

//void gf16mat_prod_16_32_asm(uint8_t *c, const uint8_t *matA, const uint8_t *b);
.global gf16mat_prod_16_32_asm
.type gf16mat_prod_16_32_asm, %function
.align 2
gf16mat_prod_16_32_asm:
    // r0 = accu
    // r1 = A
    // r2 = b
    push {r4-r11, r14}
    mov r3, #0
    .set k, 0
    .rept 4
    str r3, [r0, k]
    .set k, k+4
    .endr

    one         .req s2
    ctr1        .req s3
    ctr2        .req s4

    vmov one, #1
    vmov ctr1, #4
    1:
    ldr r3, [r2], #4      // b0 ... b31
    vmov ctr2, #8
    2:

    sbfx r11, r3, #3, #1  // b3
    sbfx r10, r3, #2, #1  // b2
    sbfx r9, r3,  #1, #1  // b1
    sbfx r8, r3,  #0, #1  // b0
    .rept 2
    ldmia r1!, {r4-r5}
    ldm r0, {r6-r7}
    accmulu32 r4, r6, r8, r9, r10, r11, r12, r14
    accmulu32 r5, r7, r8, r9, r10, r11, r12, r14
    stmia r0!, {r6-r7}
    .endr

    sub r0, r0, #16
    lsr r3, r3,  #4
    vsub.f32 ctr2, ctr2, one
    vcmp.f32 ctr2, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 2b

    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    pop {r4-r11, r14}
    bx lr
    .unreq one
    .unreq ctr1
    .unreq ctr2

.macro madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
    tst.w \b_32, #1
    itttt ne
    eorne.w \accu0, \accu0, \mat0      // out[0] ^= (b[0] & a[0])
    eorne.w \accu1, \accu1, \mat1      // out[1] ^= (b[0] & a[1])
    eorne.w \accu2, \accu2, \mat2      // out[2] ^= (b[0] & a[2])
    eorne.w \accu3, \accu3, \mat3      // out[3] ^= (b[0] & a[3])

    eor.w \tmp0, \mat0, \mat1          // tmp0 = a[0] ^ a[1]
    eor.w \tmp1, \mat2, \mat3          // tmp1 = a[2] ^ a[3]
    tst.w \b_32, #2
    itttt ne
    eorne.w \accu0, \accu0, \mat1      // out[0] ^= (b[1] & a[1])
    eorne.w \accu1, \accu1, \tmp0      // out[1] ^= (b[1] & (a[0] ^ a[1]))
    eorne.w \accu2, \accu2, \mat3      // out[2] ^= (b[1] & a[3])
    eorne.w \accu3, \accu3, \tmp1      // out[3] ^= (b[1] & (a[2] ^ a[3]))


    mov.w \tmp2, #0
    ands \tmp3, \tmp2, \b_32, LSR #3

    itttt cs
    eorcs.w \tmp2, \mat2               // tmp2 = (b[2] & a[2])
    eorcs.w \tmp3, \mat3               // tmp3 = (b[2] & a[3])
    eorcs.w \accu2, \accu2, \mat0      // out[2] ^= (b[2] & a[0])
    eorcs.w \accu3, \accu3, \mat1      // out[3] ^= (b[2] & a[1])

    tst.w \b_32, #8
    itttt ne
    eorne.w \accu2, \accu2, \mat1      // out[2] ^= (b[3] & a[1])
    eorne.w \accu3, \accu3, \tmp0      // out[3] ^= (b[3] & (a[0] ^ a[1]))
    eorne.w \tmp2, \tmp2, \mat3        // tmp2 = (b[2] & a[2]) ^ (b[3] & a[3]))
    eorne.w \tmp3, \tmp3, \tmp1        // tmp3 = (b[2] & a[3]) ^ (b[3] & (a[2] ^ a[3]))

    eor.w \accu0, \accu0, \tmp3        // out[0] ^= (b[2] & a[3]) ^ (b[3] & (a[2] ^ a[3]))
    eor.w \accu1, \accu1, \tmp2        // out[1] ^= (b[2] & a[2]) ^ (b[3] & a[3]))
    eor.w \accu1, \accu1, \tmp3        // out[1] ^= (b[2] & a[3]) ^ (b[3] & (a[2] ^ a[3]))
    eor.w \accu2, \accu2, \tmp2        // out[2] ^= (b[2] & a[2]) ^ (b[3] & a[3]))
    eor.w \accu3, \accu3, \tmp3        // out[3] ^= (b[2] & a[3]) ^ (b[3] & (a[2] ^ a[3]))
.endm


.macro bitslice out0, out1, out2, out3, in0, in1, in2, in3
    // use out3 as tmp
    and.w \out0, \in0, #0x11111111
    and.w \out3, \in1, #0x11111111
    orr.w \out0, \out0, \out3, lsl#1
    and.w \out3, \in2,  #0x11111111
    orr.w \out0, \out0, \out3, lsl#2
    and.w \out3, \in3, #0x11111111
    orr.w \out0, \out0, \out3, lsl#3

    and.w \out1, \in1, #0x22222222
    and.w \out3, \in0, #0x22222222
    orr.w \out1, \out1, \out3, lsr#1
    and.w \out3, \in2, #0x22222222
    orr.w \out1, \out1, \out3, lsl#1
    and.w \out3, \in3, #0x22222222
    orr.w \out1, \out1, \out3, lsl#2

    and.w \out2, \in2, #0x44444444
    and.w \out3, \in0, #0x44444444
    orr.w \out2, \out2, \out3, lsr#2
    and.w \out3, \in1, #0x44444444
    orr.w \out2, \out2, \out3, lsr#1
    and.w \out3, \in3, #0x44444444
    orr.w \out2, \out2, \out3, lsl#1

    and.w \out3, \in3, #0x88888888
    // in3 no longer needed; use as tmp
    and.w \in3, \in0, #0x88888888
    orr.w \out3, \out3, \in3, lsr#3
    and.w \in3, \in1, #0x88888888
    orr.w \out3, \out3, \in3, lsr#2
    and.w \in3, \in2, #0x88888888
    orr.w \out3, \out3, \in3, lsr#1
.endm

.macro bitslice_single out0, out1, out2, out3, in0
    and.w \out0, \in0, #0x11111111
    and.w \out1, \in0, #0x22222222
    lsr.w \out1, \out1, #1
    and.w \out2, \in0, #0x44444444
    lsr.w \out2, \out2, #2
    and.w \out3, \in0, #0x88888888
    lsr.w \out3, \out3, #3
.endm


.macro unbitslice_single out0, in0, in1, in2, in3
    # TODO: theoretically the masking is not needed as those elements should only have those bits set anyway
    and.w \out0, \in0, #0x11111111
    and.w \in1, \in1, #0x11111111
    orr.w \out0, \out0, \in1, lsl#1
    and.w \in2, \in2, #0x111111111
    orr.w \out0, \out0, \in2, lsl#2
    and.w \in3, \in3, #0x11111111
    orr.w \out0, \out0, \in3, lsl#3
.endm

.macro unbitslice out0, out1, out2, out3, in0, in1, in2, in3
    bitslice \out0, \out1, \out2, \out3 , \in0, \in1, \in2, \in3
.endm



//extern void gf16v_bitslice_asm(uint32_t *out, uint32_t *in, size_t len);
.global gf16v_bitslice_asm
.type gf16v_bitslice_asm, %function
.align 2
gf16v_bitslice_asm:
    push.w {r4-r9}
    1:
        ldr.w r4, [r1, #4]
        ldr.w r5, [r1, #8]
        ldr.w r6, [r1, #12]
        ldr.w r3, [r1], #16

        bitslice r7, r8, r9, r12, r3, r4, r5, r6

        str.w r8, [r0, #4]
        str.w r9, [r0, #8]
        str.w r12,[r0, #12]
        str.w r7, [r0], #16
    subs r2, r2, #1
    bne 1b
    pop.w {r4-r9}
    bx lr


//extern void gf16mat_prod_16_32_bitsliced_inner_asm(uint32_t *c, uint32_t *a, uint8_t *b);
.global gf16mat_prod_16_32_bitsliced_inner_asm
.type gf16mat_prod_16_32_bitsliced_inner_asm, %function
.align 2
gf16mat_prod_16_32_bitsliced_inner_asm:
    push {r4-r11, r14}
    c_ptr .req r0
    a_ptr .req r1
    b_ptr .req r2
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    b_32     .req r5
    tmp0    .req r6
    tmp1    .req r7
    tmp2    .req r8
    tmp3    .req r9
    accu0    .req r3
    accu1    .req r11
    accu2    .req r4
    accu3    .req r14
    one      .req s0
    ctr1     .req s1

    push.w {c_ptr}
    push.w {b_ptr}
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    vmov one, #0.5
    vmov ctr1, #2
    1:
        pop.w {b_ptr}
        ldr.w b_32, [b_ptr], #4
        push.w {b_ptr}
        .set kk, 0
        .rept 8
            ldr.w mat1, [a_ptr, #4]
            ldr.w mat2, [a_ptr, #8]
            ldr.w mat3, [a_ptr, #12]
            ldr.w mat0, [a_ptr], #16

            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b
    pop.w {b_ptr}
    pop.w {c_ptr}
    str.w accu0, [c_ptr]
    str.w accu1, [c_ptr, #4]
    str.w accu2, [c_ptr, #8]
    str.w accu3, [c_ptr, #12]

    pop.w {r4-r11, r14}
    bx lr
    .unreq c_ptr
    .unreq a_ptr
    .unreq b_ptr
    .unreq mat0
    .unreq mat1
    .unreq mat2
    .unreq mat3
    .unreq b_32
    .unreq tmp0
    .unreq tmp1
    .unreq tmp2
    .unreq tmp3
    .unreq accu0
    .unreq accu1
    .unreq accu2
    .unreq accu3
    .unreq one
    .unreq ctr1


//extern void gf16mat_prod_16_32_ontheflybitsliced_asm(uint32_t *c, uint32_t *a, uint32_t *b);
.global gf16mat_prod_16_32_ontheflybitsliced_asm
.type gf16mat_prod_16_32_ontheflybitsliced_asm, %function
.align 2
gf16mat_prod_16_32_ontheflybitsliced_asm:
    push {r4-r11, r14}
    c_ptr .req r0
    a_ptr .req r1
    b_ptr .req r2
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    b_32     .req r5
    tmp0    .req r6
    tmp1    .req r7
    tmp2    .req r8
    tmp3    .req r9
    accu0    .req r3
    accu1    .req r11
    accu2    .req r4
    accu3    .req r14

    one      .req s0
    ctr1     .req s1

    push.w {c_ptr}
    push.w {b_ptr}
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    vmov one, #0.5
    vmov ctr1, #2
    1:
        pop.w {b_ptr}
        ldr.w b_32, [b_ptr], #4
        push.w {b_ptr}
        .set kk, 0
        .rept 8
            ldr.w tmp1, [a_ptr, #4]
            ldr.w tmp2, [a_ptr, #8]
            ldr.w tmp3, [a_ptr, #12]
            ldr.w tmp0, [a_ptr], #16
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3
            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    pop.w {b_ptr}
    pop.w {c_ptr}
    // un-bitslice on the fly
    bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3
    str.w tmp0, [c_ptr]
    str.w tmp1, [c_ptr, #4]
    str.w tmp2, [c_ptr, #8]
    str.w tmp3, [c_ptr, #12]


    pop.w {r4-r11, r14}
    bx lr
    .unreq c_ptr
    .unreq a_ptr
    .unreq b_ptr
    .unreq mat0
    .unreq mat1
    .unreq mat2
    .unreq mat3
    .unreq b_32
    .unreq tmp0
    .unreq tmp1
    .unreq tmp2
    .unreq tmp3
    .unreq accu0
    .unreq accu1
    .unreq accu2
    .unreq accu3
    .unreq one
    .unreq ctr1

//extern void gf16mat_prod_16_32_ontheflybitsliced_asm(uint32_t *c, uint32_t *a, uint32_t *b);
.global gf16mat_prod_18_32_ontheflybitsliced_asm
.type gf16mat_prod_18_32_ontheflybitsliced_asm, %function
.align 2
gf16mat_prod_18_32_ontheflybitsliced_asm:
    push {r4-r11, r14}
    c_ptr .req r0
    a_ptr .req r1
    b_ptr .req r2
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    b_32     .req r5
    tmp0    .req r6
    tmp1    .req r7
    tmp2    .req r8
    tmp3    .req r9
    accu0    .req r3
    accu1    .req r11
    accu2    .req r4
    accu3    .req r14

    one      .req s0
    ctr1     .req s1

    push.w {c_ptr}
    push.w {b_ptr}
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    vmov one, #0.5
    vmov ctr1, #2
    1:
        pop.w {b_ptr}
        ldr.w b_32, [b_ptr], #4
        push.w {b_ptr}
        .set kk, 0
        .rept 8
            ldr.w tmp1, [a_ptr, #4]
            ldr.w tmp2, [a_ptr, #8]
            ldr.w tmp3, [a_ptr, #12]
            ldr.w tmp0, [a_ptr], #18
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3
            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    pop.w {b_ptr}
    pop.w {c_ptr}
    // un-bitslice on the fly
    bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3
    str.w tmp0, [c_ptr]
    str.w tmp1, [c_ptr, #4]
    str.w tmp2, [c_ptr, #8]
    str.w tmp3, [c_ptr, #12]

    // TODO: this can probably be done more efficiently.
    sub b_ptr, #16
    sub a_ptr, #8*4*18-16
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    push.w {c_ptr}
    push.w {b_ptr}

    vmov one, #0.5
    vmov ctr1, #2
    1:
        pop.w {b_ptr}
        ldr.w b_32, [b_ptr], #4
        push.w {b_ptr}
        .set kk, 0
        .rept 8
            ldrh.w tmp0, [a_ptr], #18
            // bitslice on the fly
            bitslice_single mat0, mat1, mat2, mat3, tmp0
            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    pop.w {b_ptr}
    pop.w {c_ptr}
    unbitslice_single tmp0, accu0, accu1, accu2, accu3
    strh.w tmp0, [c_ptr, #16]

    pop.w {r4-r11, r14}
    bx lr
    .unreq c_ptr
    .unreq a_ptr
    .unreq b_ptr
    .unreq mat0
    .unreq mat1
    .unreq mat2
    .unreq mat3
    .unreq b_32
    .unreq tmp0
    .unreq tmp1
    .unreq tmp2
    .unreq tmp3
    .unreq accu0
    .unreq accu1
    .unreq accu2
    .unreq accu3
    .unreq one
    .unreq ctr1

//extern void gf16mat_prod_512_32_bitsliced_inner_asm(uint32_t *c, uint32_t *a, uint8_t *b);
.global gf16mat_prod_512_32_bitsliced_inner_asm
.type gf16mat_prod_512_32_bitsliced_inner_asm, %function
.align 2
gf16mat_prod_512_32_bitsliced_inner_asm:
    push {r4-r11, r14}
    c_ptr .req r0
    a_ptr .req r1
    b_ptr .req r2
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    b_32     .req r5
    tmp0    .req r6
    tmp1    .req r7
    tmp2    .req r8
    tmp3    .req r9
    accu0    .req r3
    accu1    .req r11
    accu2    .req r4
    accu3    .req r14

    one      .req s0
    ctr1     .req s1
    ctr2     .req s2

    push.w {c_ptr}
    push.w {b_ptr}

    vmov one, #0.5
    vmov ctr1, #16
    1:
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    vmov ctr2, #2
    2:
        //.if ii != 0
        pop.w {b_ptr}
        //.endif
        ldr.w b_32, [b_ptr], #4
        //.if ii != 3
        push.w {b_ptr}
        //.endif
        .set kk, 0
        .rept 8
            ldr.w mat0, [a_ptr]
            ldr.w mat1, [a_ptr, #4]
            ldr.w mat2, [a_ptr, #8]
            ldr.w mat3, [a_ptr, #12]
            add a_ptr, #512

            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
        .endr
    //.endr
    vsub.f32 ctr2, ctr2, one
    vcmp.f32 ctr2, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 2b

    pop.w {b_ptr}
    pop.w {c_ptr}
    str.w accu1, [c_ptr, #4]
    str.w accu2, [c_ptr, #8]
    str.w accu3, [c_ptr, #12]
    str.w accu0, [c_ptr], #16
    push {c_ptr}
    mov r3, #16368 //(32*512-16)
    sub a_ptr, a_ptr, r3
    sub b_ptr, b_ptr, #16
    push.w {b_ptr}
    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    add sp, #8
    pop.w {r4-r11, r14}
    bx lr
    .unreq c_ptr
    .unreq a_ptr
    .unreq b_ptr
    .unreq mat0
    .unreq mat1
    .unreq mat2
    .unreq mat3
    .unreq b_32
    .unreq tmp0
    .unreq tmp1
    .unreq tmp2
    .unreq tmp3
    .unreq accu0
    .unreq accu1
    .unreq accu2
    .unreq accu3

//extern void gf16mat_prod_512_32_ontheflybitsliced_asm(uint32_t *c, uint32_t *a, uint8_t *b);
.global gf16mat_prod_512_32_ontheflybitsliced_asm
.type gf16mat_prod_512_32_bitsliced_inner_asm, %function
.align 2
gf16mat_prod_512_32_ontheflybitsliced_asm:
    push {r4-r11, r14}
    c_ptr .req r0
    a_ptr .req r1
    b_ptr .req r2
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    b_32     .req r5
    tmp0    .req r6
    tmp1    .req r7
    tmp2    .req r8
    tmp3    .req r9
    accu0    .req r3
    accu1    .req r11
    accu2    .req r4
    accu3    .req r14


    one      .req s0
    ctr1     .req s1
    ctr2     .req s2

    push.w {c_ptr}
    push.w {b_ptr}

    vmov one, #0.5
    vmov ctr1, #16
    1:
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    vmov ctr2, #2
    2:
        pop.w {b_ptr}
        ldr.w b_32, [b_ptr], #4
        push.w {b_ptr}
        .set kk, 0
        .rept 8
            ldr.w tmp0, [a_ptr]
            ldr.w tmp1, [a_ptr, #4]
            ldr.w tmp2, [a_ptr, #8]
            ldr.w tmp3, [a_ptr, #12]
            add a_ptr, #512
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3

            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr2, ctr2, one
    vcmp.f32 ctr2, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 2b


    pop.w {b_ptr}
    pop.w {c_ptr}
    // un-bitslice on the fly
    bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3

    str.w tmp1, [c_ptr, #4]
    str.w tmp2, [c_ptr, #8]
    str.w tmp3, [c_ptr, #12]
    str.w tmp0, [c_ptr], #16
    push {c_ptr}
    mov r3, #16368 //(32*512-16)
    sub a_ptr, a_ptr, r3
    sub b_ptr, b_ptr, #16
    push.w {b_ptr}

    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    add sp, #8

    pop.w {r4-r11, r14}
    bx lr
    .unreq c_ptr
    .unreq a_ptr
    .unreq b_ptr
    .unreq mat0
    .unreq mat1
    .unreq mat2
    .unreq mat3
    .unreq b_32
    .unreq tmp0
    .unreq tmp1
    .unreq tmp2
    .unreq tmp3
    .unreq accu0
    .unreq accu1
    .unreq accu2
    .unreq accu3
    .unreq one
    .unreq ctr1
    .unreq ctr2

//batch_quad_trimat_eval_gf16_32_16( unsigned char * y, const unsigned char * trimat, const unsigned char * x,);
//.global batch_quad_trimat_eval_gf16_32_16
//.type batch_quad_trimat_eval_gf16_32_16, %function
//.align 2
//batch_quad_trimat_eval_gf16_32_16:
//    push {r4-r11, r14}
//    y      .req r0  // batch_size bytes -> 16 bytes
 //   x      .req r1  // dim/2 bytes -> 16 bytes
 //   trimat .req r2  // dim*dim*batch_size bytes -> 16384 bytes

 //   tmp0   .req r3
 //   tmp1   .req r4
 //   tmp2   .req r5
 //   tmp3   .req r6

    // dim = 32
    // batch_size = 16
 //   .set dim, 32

    /* stack:
            _x  -> dim bytes
     */
//    sub sp, #dim
//    .set k, 0
//    .rept dim/4
//        ldr r3, [x, #k*4]
//        strb r3, [sp, #k*8]
 //       .set j, 1
 //       .rept 7
 //           lsr r3, r3, #4
 //           strb r3, [sp, #k*8+j]
 //           .set j, j+1
///       .endr
 //       .set k, k+1
  //  .endr




//    add sp, #dim
//    pop {r4-r11, r14}
//    bx lr

.macro accmulu32_new a, c, x, t1, t2, t3

    sbfx \t1, \x,  #0, #1 //b0

    and \t2, \a, \t1  // a&b0
    eor \c, \c, \t2
    sbfx \t3, \x,  #2, #1 //b2
    and \t1, \a, \t3  // a&b2
    mul2 \a, \t2

    sbfx \t3, \x,  #1, #1 //b1
    and \t2, \a, \t3 // 2a & b1
    sbfx \t3, \x,  #3, #1 //b3
    and \a, \a, \t3  // 2a & b3
    eor \c, \c, \t2
    eor \a, \t1, \a

    lsl \t1, \a, #2  // axb1
    and \t1, \t1, #0xcccccccc //a0b1
    eor \c, \c, \t1
    and \a, \a, #0xcccccccc //a1b1
    eor \c, \c, \a

    lsr \a, \a, #2 // a1b1_2

    //gf4v_mul_2_u32
    mul2 \a, \t1
    eor \c, \c, \a
.endm

//void batch_quad_trimat_eval_gf16_32_16_inner_asm(unsigned char *tmp, unsigned char *trimat, unsigned char *x, unsigned int i)
.global batch_quad_trimat_eval_gf16_32_16_inner_asm
.type batch_quad_trimat_eval_gf16_32_16_inner_asm, %function
.align 2
batch_quad_trimat_eval_gf16_32_16_inner_asm:
    push {r4-r11, r14}

    tmp     .req r0
    trimat  .req r1
    x_ptr   .req r2
    x       .req r4
    b0      .req r4
    b1      .req r5
    b2      .req r6
    b3      .req r7
    accu0    .req r8
    accu1    .req r9
    accu2    .req r10
    accu3    .req r11
    a0      .req r12
    t1      .req r14
    t2      .req r0
    t3      .req r5


    add x_ptr, x_ptr, r3
    mov r4, #32
    sub r3, r4, r3

    push {tmp}
    //TODO: inline this
    mov accu0, #0
    mov accu1, #0
    mov accu2, #0
    mov accu3, #0

    1:
    // TODO load full word and then split it up
    ldrb x, [x_ptr], #1
    //sbfx b3, x,  #3, #1
    //sbfx b2, x,  #2, #1
    //sbfx b1, x,  #1, #1
    //sbfx b0, x,  #0, #1
    ldr a0, [trimat], #4
    accmulu32_new a0, accu0, x, t1, t2, t3
    ldr a0, [trimat], #4
    accmulu32_new a0, accu1, x, t1, t2, t3
    ldr a0, [trimat], #4
    accmulu32_new a0, accu2, x, t1, t2, t3
    ldr a0, [trimat], #4
    accmulu32_new a0, accu3, x, t1, t2, t3
    subs r3, r3, #1
    bne 1b

    pop {tmp}
    str accu0, [tmp, #0]
    str accu1, [tmp, #4]
    str accu2, [tmp, #8]
    str accu3, [tmp, #12]

    pop {r4-r11, r14}
    bx lr


    .unreq tmp
    .unreq trimat
    .unreq x_ptr
    .unreq x
    .unreq b0
    .unreq b1
    .unreq b2
    .unreq b3
    .unreq accu0
    .unreq accu1
    .unreq accu2
    .unreq accu3
    .unreq a0
    .unreq t1
    .unreq t2


//void batch_quad_trimat_eval_gf16_32_16_inner2_asm(unsigned *y, unsigned char *trimat, unsigned char *_x)
.global batch_quad_trimat_eval_gf16_32_16_inner2_asm
.type batch_quad_trimat_eval_gf16_32_16_inner2_asm, %function
.align 2
batch_quad_trimat_eval_gf16_32_16_inner2_asm:
    push {r4-r11, r14}
    y       .req r4
    tmp     .req r0
    trimat  .req r1
    x_ptr   .req r2
    ctr     .req r3
    y0      .req r5
    t1      .req r6
    t2      .req r7
    tmp0      .req r8

    b0      .req r9
    b1      .req r10
    b2      .req r11
    b3      .req r12
    x       .req r9

    mov y, r0

    //gf256v_set_zero( y , size_batch);
    //TODO: inline this
    mov r12, #0
    str r12, [y]
    str r12, [y,#4]
    str r12, [y,#8]
    str r12, [y,#12]

    //unsigned char tmp[16];
    sub sp, sp, #16
    mov tmp, sp

    mov ctr, #0
    1:
        push {tmp, trimat, x_ptr, ctr}
        bl batch_quad_trimat_eval_gf16_32_16_inner_asm
        pop {tmp, trimat, x_ptr, ctr}
        mov r5, #32
        sub r5, r5, ctr
        add trimat, trimat, r5, lsl#4

        ldrb x, [x_ptr, ctr]
        sbfx b3, x,  #3, #1
        sbfx b2, x,  #2, #1
        sbfx b1, x,  #1, #1
        sbfx b0, x,  #0, #1

        ldr tmp0, [tmp, #0]
        ldr y0, [y]
        accmulu32 tmp, y0, b0, b1, b2, b3, t1, t2
        str y0, [y]

        ldr tmp0, [tmp, #4]
        ldr y0, [y, #4]
        accmulu32 tmp, y0, b0, b1, b2, b3, t1, t2
        str y0, [y, #4]

        ldr tmp0, [tmp, #8]
        ldr y0, [y, #8]
        accmulu32 tmp, y0, b0, b1, b2, b3, t1, t2
        str y0, [y, #8]

        ldr tmp0, [tmp, #12]
        ldr y0, [y, #12]
        accmulu32 tmp, y0, b0, b1, b2, b3, t1, t2
        str y0, [y, #12]

    add ctr, ctr, #1
    cmp ctr, #32
    bne 1b

    add sp, #16
    pop {r4-r11, r14}
    bx lr

.unreq tmp0

//extern void gf16mat_prod_512_32_ontheflybitsliced_asm(uint32_t *c, uint32_t *a, uint8_t *b);
.global gf16mat_prod_512_36_ontheflybitsliced_asm
.type gf16mat_prod_512_36_bitsliced_inner_asm, %function
.align 2
gf16mat_prod_512_36_ontheflybitsliced_asm:
    push {r4-r11, r14}
    c_ptr .req r0
    a_ptr .req r1
    b_ptr .req r2
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    b_32     .req r5
    tmp0    .req r6
    tmp1    .req r7
    tmp2    .req r8
    tmp3    .req r9
    accu0    .req r3
    accu1    .req r11
    accu2    .req r4
    accu3    .req r14


    one      .req s0
    ctr1     .req s1
    ctr2     .req s2

    push.w {c_ptr}
    push.w {b_ptr}

    vmov one, #0.5
    vmov ctr1, #16
    1:
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    vmov ctr2, #2
    2:
        pop.w {b_ptr}
        ldr.w b_32, [b_ptr], #4
        push.w {b_ptr}
        .set kk, 0
        .rept 8
            ldr.w tmp0, [a_ptr]
            ldr.w tmp1, [a_ptr, #4]
            ldr.w tmp2, [a_ptr, #8]
            ldr.w tmp3, [a_ptr, #12]
            add a_ptr, #512
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3

            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr2, ctr2, one
    vcmp.f32 ctr2, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 2b

    pop.w {b_ptr}
    ldrh.w b_32, [b_ptr]
    push.w {b_ptr}
    .set kk, 0
    .rept 4
        ldr.w tmp0, [a_ptr]
        ldr.w tmp1, [a_ptr, #4]
        ldr.w tmp2, [a_ptr, #8]
        ldr.w tmp3, [a_ptr, #12]
        //TODO: do this smarter
        add a_ptr, #512
        // bitslice on the fly
        bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3
        madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
        .if kk != 3
            lsr.w b_32, b_32, #4
        .endif
        .set kk, kk+1
    .endr
    pop.w {b_ptr}
    pop.w {c_ptr}
    // un-bitslice on the fly
    bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3

    str.w tmp1, [c_ptr, #4]
    str.w tmp2, [c_ptr, #8]
    str.w tmp3, [c_ptr, #12]
    str.w tmp0, [c_ptr], #16
    push {c_ptr}
    mov r3, #18416 // 36*512-16)
    sub a_ptr, a_ptr, r3
    sub b_ptr, b_ptr, #16
    push.w {b_ptr}

    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    add sp, #8

    pop.w {r4-r11, r14}
    bx lr
    .unreq c_ptr
    .unreq a_ptr
    .unreq b_ptr
    .unreq mat0
    .unreq mat1
    .unreq mat2
    .unreq mat3
    .unreq b_32
    .unreq tmp0
    .unreq tmp1
    .unreq tmp2
    .unreq tmp3
    .unreq accu0
    .unreq accu1
    .unreq accu2
    .unreq accu3
    .unreq one
    .unreq ctr1
    .unreq ctr2
