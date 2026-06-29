// SPDX-License-Identifier: MIT
// NTRU+864 (X^3-zeta, degree-3) Cortex-M4 packed basemul (q=3457). Pairings
// host-verified in results/ntruplus_opt/basemul_check3.c (bit-exact vs scalar C
// over 5M inputs). Packs a/b as A01=a0|a1, A12=a1|a2, B01=b0|b1, B12=b1|b2.
//   void basemul(int16_t r[3], const int16_t a[3], const int16_t b[3], int16_t zeta);
// QINV=12929 Q=3457 RSQ=867 R=-147
//
// h0 = a1*b2 + a2*b1 = smuadx(A12,B12);  h1 = a2*b2 = smultt(A12,B12)
// r0 = mont(h0r*zeta + a0*b0)
// r1 = mont(h1r*zeta + a0*b1 + a1*b0)
// r2 = mont(a2*b0 + a1*b1 + a0*b2)
// r[i] = mont(r_i * RSQ)   [basemul]   or   mont(c[i]*R + r_i*RSQ)  [basemul_add]

.syntax unified
.cpu cortex-m4
.thumb

// montgomery_reduce of int32 \x in place; r8=qinv, r9=q, clobbers r2.
.macro MONT x
    mul    r2, \x, r8
    sxth   r2, r2
    smulbb r2, r2, r9
    sub    \x, \x, r2
    asr    \x, \x, #16
.endm

.global basemul
.type   basemul, %function
.align 2
basemul:
    push   {r4-r11, lr}
    ldr    r4, [r1]                 // A01 = a0|a1
    ldr    r5, [r1, #2]             // A12 = a1|a2  (unaligned)
    ldr    r6, [r2]                 // B01 = b0|b1
    ldr    r7, [r2, #2]             // B12 = b1|b2  (unaligned)
    movw   r8, #12929               // qinv
    movw   r9, #3457                // q
    // stage 1: high-conv sums -> mont
    smuadx r10, r5, r7              // h0 = a1*b2 + a2*b1
    MONT   r10                      // h0r
    smultt r11, r5, r7              // h1 = a2*b2
    MONT   r11                      // h1r
    // stage 2
    smulbb r12, r10, r3             // h0r*zeta
    smlabb r12, r4, r6, r12         // + a0*b0
    MONT   r12                      // out r0
    smuadx lr, r4, r6               // a0*b1 + a1*b0
    smlabb lr, r11, r3, lr          // + h1r*zeta
    MONT   lr                       // out r1
    smultb r11, r5, r6              // a2*b0
    smlatt r11, r4, r6, r11         // + a1*b1
    smlabt r11, r4, r7, r11         // + a0*b2
    MONT   r11                      // out r2
    // [FIX-ASM] clean convention: basemul output is R^-1 domain (a*b*R^-1).
    // Removed trailing mont(out_ri * RSQ) that lifted result by extra factor R
    // (asm = clean*R). Store out_ri directly -> bytes match clean.
    strh   r12, [r0]
    strh   lr,  [r0, #2]
    strh   r11, [r0, #4]
    pop    {r4-r11, pc}
.size basemul, .-basemul

// void basemul_add(int16_t r[3], a[3], b[3], const int16_t c[3], int16_t zeta);
// same stage1/2 as basemul; final r[i] = mont(c[i]*R + out_ri*RSQ).  R=-147
.global basemul_add
.type   basemul_add, %function
.align 2
basemul_add:
    push   {r4-r11, lr}
    ldr    r4, [r1]                 // A01
    ldr    r5, [r1, #2]             // A12
    ldr    r6, [r2]                 // B01
    ldr    r7, [r2, #2]             // B12
    ldr    r1, [sp, #36]            // zeta (5th arg) -> r1
    movw   r8, #12929               // qinv
    movw   r9, #3457                // q
    smuadx r10, r5, r7              // h0 = a1*b2 + a2*b1
    MONT   r10
    smultt r11, r5, r7              // h1 = a2*b2
    MONT   r11
    // [FIX-ASM] clean convention: basemul_add = mont(c[i]*(-147) + raw stage-2
    // accumulator) = c[i] + a*b*R^-1. Previous code did mont(c*R + out_ri*RSQ)
    // = c + out_ri*R (extra factor R). Fix: fold c[i]*R into each lane's
    // accumulator BEFORE its single MONT. R=-147 (u16 65389) kept on stack since
    // r4..r7 hold a/b packs; c pointer in r3.
    movw   r2, #65389               // R = -147 (low16)
    push   {r2}                     // [sp] = R
    smulbb r12, r10, r1             // h0r*zeta
    smlabb r12, r4, r6, r12         // + a0*b0
    ldrsh  r2, [r3]                 // c0
    ldr    r10, [sp]                // R  (h0r dead)
    smlabb r12, r2, r10, r12        // + c0*R
    MONT   r12                      // r[0]
    strh   r12, [r0]
    smuadx lr, r4, r6               // a0*b1 + a1*b0
    smlabb lr, r11, r1, lr          // + h1r*zeta
    ldrsh  r2, [r3, #2]            // c1
    ldr    r10, [sp]                // R
    smlabb lr, r2, r10, lr          // + c1*R
    MONT   lr                       // r[1]
    strh   lr,  [r0, #2]
    smultb r11, r5, r6              // a2*b0
    smlatt r11, r4, r6, r11         // + a1*b1
    smlabt r11, r4, r7, r11         // + a0*b2
    ldrsh  r2, [r3, #4]            // c2
    ldr    r10, [sp]                // R
    smlabb r11, r2, r10, r11        // + c2*R
    MONT   r11                      // r[2]
    strh   r11, [r0, #4]
    add    sp, sp, #4               // pop reserved R slot
    pop    {r4-r11, pc}
.size basemul_add, .-basemul_add
