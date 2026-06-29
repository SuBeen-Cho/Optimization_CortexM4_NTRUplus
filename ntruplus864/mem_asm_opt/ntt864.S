// SPDX-License-Identifier: MIT
// NTRU+864 Cortex-M4: radix-2 Cooley-Tukey layers of ntt() (q=3457).
// Unlike 768, the 1152/864 radix-2 layers apply barrett_reduce to BOTH outputs
// inside the butterfly (no separate final pass). steps 24,12,6,3 (base X^3).
//
// PACKED-16: two adjacent butterflies (i,i+1) share zeta within a block. Load
// r[i]|r[i+1] and r[i+step]|r[i+step+1] as words, Montgomery-multiply both halves
// of the second by zeta (SMULBB/SMULBT), repack with PKHBT, do BOTH add/sub with
// SADD16/SSUB16, then barrett each lane and repack. An odd-step scalar tail (used
// for N=864's step=3) handles the lone butterfly per block. Bit-identical, KAT-gated.
//   void ntt_ct_layers(int16_t *r, const int16_t *zp);  zp=&zetas[18]
//
// montgomery_reduce(p): m=(int16)(p*QINV); (p - m*Q) >> 16     QINV=12929 Q=3457
// barrett_reduce(a):     a - (((v*a + (1<<25)) >> 26) * Q)     v=19412

.syntax unified
.cpu cortex-m4
.thumb

.equ NCOEFF,     864
.equ FIRST_STEP, 24
.equ LAST_STEP,  3

.global ntt_ct_layers
.type   ntt_ct_layers, %function
.align 2
ntt_ct_layers:
    push   {r4-r11, lr}
    movw   r2, #12929              // qinv
    movw   r3, #3457               // q
    mov    r4, #FIRST_STEP         // step (coeffs)
.Llayer:
    lsl    r5, r4, #1              // step_bytes = step*2
    mov    r6, #0                  // start_byte
.Lstart:
    ldrsh  r7, [r1], #2            // zeta = *zp++
    add    r8, r0, r6              // pA = &r[start]
    add    r9, r8, r5              // pB = pA + step_bytes
    asr    r10, r4, #1             // cnt = step/2 (packed butterflies)
    push   {r5, r6}                // free r5,r6 as butterfly scratch
.Lbfly:
    ldr    r5, [r8]                // A = r[i]|r[i+1]
    ldr    r6, [r9]                // B = r[i+step]|r[i+step+1]
    smulbt r11, r7, r6             // zeta*B.hi
    smulbb r6,  r7, r6             // zeta*B.lo
    mul    r12, r6, r2             // mont(B.lo prod)
    sxth   r12, r12
    smulbb r12, r12, r3
    sub    r6,  r6, r12
    asr    r6,  r6, #16            // t_lo
    mul    r12, r11, r2            // mont(B.hi prod)
    sxth   r12, r12
    smulbb r12, r12, r3
    sub    r11, r11, r12
    asr    r11, r11, #16           // t_hi
    pkhbt  r6, r6, r11, lsl #16    // T = t_lo|t_hi
    sadd16 r11, r5, r6             // sum  = A + T   -> r[i]|r[i+1]
    ssub16 r5,  r5, r6             // diff = A - T   -> r[i+step]|..
    // barrett(sum) -> store to pA
    movw   lr, #19412              // v
    smulbb r12, lr, r11            // v*sum.lo
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    sxth   r6, r11                 // sum.lo signed
    sub    r6, r6, r12             // barrett(sum.lo)
    smulbt r12, lr, r11            // v*sum.hi
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    asr    r11, r11, #16           // sum.hi signed
    sub    r11, r11, r12           // barrett(sum.hi)
    pkhbt  r6, r6, r11, lsl #16
    str    r6, [r8], #4            // r[i]|r[i+1]
    // barrett(diff) -> store to pB
    smulbb r12, lr, r5             // v*diff.lo
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    sxth   r6, r5                  // diff.lo signed
    sub    r6, r6, r12             // barrett(diff.lo)
    smulbt r12, lr, r5             // v*diff.hi
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    asr    r5, r5, #16             // diff.hi signed
    sub    r5, r5, r12             // barrett(diff.hi)
    pkhbt  r6, r6, r5, lsl #16
    str    r6, [r9], #4            // r[i+step]|..
    subs   r10, r10, #1
    bne    .Lbfly
    // odd-step tail (step=3 for N=864): one lone scalar barrett-inside butterfly
    tst    r4, #1
    beq    .Lnotail
    ldrsh  r6, [r9]                // cB
    smulbb r6, r7, r6              // zeta*cB
    mul    r12, r6, r2
    sxth   r12, r12
    smulbb r12, r12, r3
    sub    r6, r6, r12
    asr    r6, r6, #16             // t1
    ldrsh  r5, [r8]                // cA
    sub    r11, r5, r6             // cA - t1
    add    r5, r5, r6              // cA + t1
    movw   lr, #19412
    smulbb r12, lr, r11
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    sub    r11, r11, r12           // barrett(cA - t1)
    strh   r11, [r9]
    smulbb r12, lr, r5
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    sub    r5, r5, r12             // barrett(cA + t1)
    strh   r5, [r8]
.Lnotail:
    pop    {r5, r6}                // restore step_bytes, start_byte
    add    r6, r6, r5, lsl #1      // start_byte += 2*step_bytes
    movw   r11, #(NCOEFF*2)
    cmp    r6, r11
    blt    .Lstart
    asr    r4, r4, #1              // step >>= 1
    cmp    r4, #LAST_STEP
    bge    .Llayer
    pop    {r4-r11, pc}
.size ntt_ct_layers, .-ntt_ct_layers
