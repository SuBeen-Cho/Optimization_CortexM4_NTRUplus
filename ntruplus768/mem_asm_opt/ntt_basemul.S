// SPDX-License-Identifier: MIT
// NTRU+ (X^4-zeta) Cortex-M4 packed basemul (q=3457). Pairings host-verified in
// results/ntruplus_opt/basemul_check.c (bit-exact vs scalar C over 5M inputs).
// Packs a/b as A01=a0|a1, A23=a2|a3, B01, B23 and uses SMUAD/SMUADX/SMLA.
//   void basemul(int16_t r[4], const int16_t a[4], const int16_t b[4], int16_t zeta);
// QINV=12929 Q=3457 RSQ=867

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
    ldr    r5, [r1, #4]             // A23 = a2|a3
    ldr    r6, [r2]                 // B01 = b0|b1
    ldr    r7, [r2, #4]             // B23 = b2|b3
    movw   r8, #12929               // qinv
    movw   r9, #3457                // q
    // stage 1: high-conv sums -> mont
    smultt r10, r5, r7              // h2 = a3*b3
    MONT   r10                      // h2r
    smuadx r11, r5, r7              // h1 = a2*b3 + a3*b2
    MONT   r11                      // h1r
    smulbb r12, r5, r7              // a2*b2
    smlatt r12, r5, r6, r12         // + a3*b1
    smlatt r12, r4, r7, r12         // + a1*b3
    MONT   r12                      // h0r
    // stage 2: ri = mont(hir*zeta + low-conv)
    smulbb lr, r12, r3              // h0r*zeta
    smlabb lr, r4, r6, lr           // + a0*b0
    MONT   lr                       // out r0
    smuadx r12, r4, r6              // a0*b1 + a1*b0
    smlabb r12, r11, r3, r12        // + h1r*zeta
    MONT   r12                      // out r1
    smulbb r11, r4, r7              // a0*b2
    smlatt r11, r4, r6, r11         // + a1*b1
    smlabb r11, r5, r6, r11         // + a2*b0
    smlabb r11, r10, r3, r11        // + h2r*zeta
    MONT   r11                      // out r2
    smuadx r10, r4, r7              // a0*b3 + a1*b2
    smuadx r1,  r5, r6              // a2*b1 + a3*b0  (r1 free: a ptr done)
    add    r10, r10, r1             // x3
    MONT   r10                      // out r3
    // [FIX-ASM] clean convention: basemul output is R^-1 domain (a*b*R^-1).
    // The previous trailing mont(out_ri * RSQ) lifted the result by an extra
    // factor R (asm = clean*R), diverging from clean serialization. Removed —
    // out_ri == clean basemul output, stored directly so bytes match clean.
    strh   lr,  [r0]
    strh   r12, [r0, #2]
    strh   r11, [r0, #4]
    strh   r10, [r0, #6]
    pop    {r4-r11, pc}
.size basemul, .-basemul

// void basemul_add(int16_t r[4], a[4], b[4], const int16_t c[4], int16_t zeta);
// same stage1/2 as basemul; final r[i] = mont(c[i]*R + out_ri*RSQ).  R=-147
.global basemul_add
.type   basemul_add, %function
.align 2
basemul_add:
    push   {r4-r11, lr}
    ldr    r4, [r1]                 // A01
    ldr    r5, [r1, #4]             // A23
    ldr    r6, [r2]                 // B01
    ldr    r7, [r2, #4]             // B23
    ldr    r1, [sp, #36]            // zeta (5th arg) -> r1
    movw   r8, #12929               // qinv
    movw   r9, #3457                // q
    smultt r10, r5, r7              // h2 = a3*b3
    MONT   r10
    smuadx r11, r5, r7              // h1 = a2*b3 + a3*b2
    MONT   r11
    smulbb r12, r5, r7              // a2*b2
    smlatt r12, r5, r6, r12         // + a3*b1
    smlatt r12, r4, r7, r12         // + a1*b3
    MONT   r12                      // h0r
    // [FIX-ASM] clean convention: basemul_add computes
    //   r[i] = mont(c[i]*(-147) + <raw stage-2 accumulator>)
    // i.e. result = c[i] + a*b*R^-1 (R^-1 product domain). The previous code
    // reduced the accumulator (out_ri) first then did mont(c*R + out_ri*RSQ)
    // = c + out_ri*R, lifting the product by an extra factor R vs clean. Fix:
    // build each lane's stage-2 accumulator with the SAME a*b math (r4=A01,
    // r5=A23, r6=B01, r7=B23 stay live), then fold c[i]*R into it BEFORE its
    // single MONT and store directly. R=-147 (u16 65389) is kept on the stack
    // (one reserved slot) because r4..r7 are all occupied by a/b packs.
    movw   r2, #65389               // R = -147 (low16)
    push   {r2}                     // [sp] = R  (stack stays 8-byte aligned: 9+1 words)
    // r[0] = mont(c0*R + h0r*zeta + a0*b0)
    smulbb lr, r12, r1              // h0r*zeta
    smlabb lr, r4, r6, lr           // + a0*b0
    ldrsh  r2, [r3]                 // c0
    ldr    r12, [sp]                // R
    smlabb lr, r2, r12, lr          // + c0*R
    MONT   lr                       // r[0]
    strh   lr,  [r0]
    // r[1] = mont(c1*R + (a0*b1+a1*b0) + h1r*zeta)
    smuadx r12, r4, r6              // a0*b1 + a1*b0
    smlabb r12, r11, r1, r12        // + h1r*zeta
    ldrsh  r2, [r3, #2]            // c1
    ldr    lr, [sp]                 // R
    smlabb r12, r2, lr, r12         // + c1*R
    MONT   r12                      // r[1]
    strh   r12, [r0, #2]
    // r[2] = mont(c2*R + a0*b2 + a1*b1 + a2*b0 + h2r*zeta)
    smulbb r11, r4, r7              // a0*b2
    smlatt r11, r4, r6, r11         // + a1*b1
    smlabb r11, r5, r6, r11         // + a2*b0
    smlabb r11, r10, r1, r11        // + h2r*zeta
    ldrsh  r2, [r3, #4]            // c2
    ldr    lr, [sp]                 // R
    smlabb r11, r2, lr, r11         // + c2*R
    MONT   r11                      // r[2]
    strh   r11, [r0, #4]
    // r[3] = mont(c3*R + a0*b3 + a1*b2 + a2*b1 + a3*b0)
    smuadx r10, r4, r7              // a0*b3 + a1*b2
    smuadx r4,  r5, r6              // a2*b1 + a3*b0  (r4/A01 dead now)
    add    r10, r10, r4             // x3
    ldrsh  r2, [r3, #6]            // c3
    ldr    lr, [sp]                 // R
    smlabb r10, r2, lr, r10         // + c3*R
    MONT   r10                      // r[3]
    strh   r10, [r0, #6]
    add    sp, sp, #4               // pop reserved R slot
    pop    {r4-r11, pc}
.size basemul_add, .-basemul_add

