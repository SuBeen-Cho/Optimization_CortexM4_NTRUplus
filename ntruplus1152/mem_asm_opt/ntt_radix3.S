// SPDX-License-Identifier: MIT
// NTRU+ Cortex-M4: one radix-3 (cube-root OMEGA) NTT layer (q=3457).
// Generic over the three sets (q, OMEGA identical; only `step` differs):
//   768 step=128 (1 layer); 1152 step=192,64; 864 step=144,48.   All steps EVEN.
// For each start in [0,n) step 3*step, with zeta1,zeta2 = *zp++,*zp++:
//   for i in [start, start+step):
//     t1 = mont(zeta1, r[i+step]); t2 = mont(zeta2, r[i+2*step]);
//     t3 = mont(OMEGA, t1 - t2);
//     r[i+2*step] = r[i] - t1 - t3;  r[i+step] = r[i] - t2 + t3;  r[i] = r[i] + t1 + t2;
//
// PACKED-16: two adjacent coeffs (i,i+1) share (zeta1,zeta2,OMEGA) within a block.
// Load A=r[i]|r[i+1], B=r[i+s]|.., C=r[i+2s]|.. as words; the three Montgomery
// multiplies stay dual-scalar (SMULBB low lane / SMULBT high lane, twiddle in low
// halfword), each product reduced and PKHBT-packed; all additive combines use one
// SADD16/SSUB16 over both lanes. step is even for every layer -> step/2 packed iters,
// no tail. Bit-identical to the C / scalar (host-verified 5M in radix3_pack_check.c);
// KAT-gated.   void ntt_radix3_layer(int16_t *r, const int16_t *zp, int n, int step);
//
// montgomery_reduce(p): m=(int16)(p*QINV); (p - m*Q) >> 16     QINV=12929 Q=3457 OMEGA=-886

.syntax unified
.cpu cortex-m4
.thumb

.global ntt_radix3_layer
.type   ntt_radix3_layer, %function
.align 2
ntt_radix3_layer:
    push   {r4-r11, lr}
    sub    sp, sp, #20
    str    r0, [sp, #0]             // r_base
    lsl    r4, r3, #1               // step_bytes = step*2
    str    r4, [sp, #4]
    add    r5, r4, r4, lsl #1       // 3*step_bytes
    str    r5, [sp, #8]
    lsl    r2, r2, #1               // n_bytes = n*2
    str    r2, [sp, #12]
    mov    r6, #0
    str    r6, [sp, #16]            // start_byte = 0
    movw   r2, #12929               // qinv (persistent)
    movw   r3, #3457                // q (persistent)
.Lr3_start:
    ldr    r10, [sp, #0]            // r_base
    ldr    r11, [sp, #16]           // start_byte
    add    r10, r10, r11            // pA = &r[start]
    ldr    r4, [sp, #4]             // step_bytes
    add    r8, r10, r4              // pB = pA + step_bytes
    add    r9, r10, r4, lsl #1      // pC = pA + 2*step_bytes
    ldrsh  r6, [r1], #2             // zeta1
    ldrsh  r7, [r1], #2             // zeta2
    lsr    r5, r4, #2               // cnt = step/2  (packed iters; r4=step_bytes)
.Lr3_bfly:
    // ---- T1 = mont(zeta1 * B), B = r[i+step]|r[i+step+1] ----
    ldr    r4, [r8]                 // B
    smulbt lr, r6, r4               // z1*B.hi
    smulbb r4, r6, r4               // z1*B.lo
    mul    r0, r4, r2
    sxth   r0, r0
    smulbb r0, r0, r3
    sub    r4, r4, r0
    asr    r4, r4, #16              // t1.lo
    mul    r0, lr, r2
    sxth   r0, r0
    smulbb r0, r0, r3
    sub    lr, lr, r0
    asr    lr, lr, #16              // t1.hi
    pkhbt  r11, r4, lr, lsl #16     // T1
    // ---- T2 = mont(zeta2 * C), C = r[i+2step]|.. ----
    ldr    r4, [r9]                 // C
    smulbt lr, r7, r4               // z2*C.hi
    smulbb r4, r7, r4               // z2*C.lo
    mul    r0, r4, r2
    sxth   r0, r0
    smulbb r0, r0, r3
    sub    r4, r4, r0
    asr    r4, r4, #16              // t2.lo
    mul    r0, lr, r2
    sxth   r0, r0
    smulbb r0, r0, r3
    sub    lr, lr, r0
    asr    lr, lr, #16              // t2.hi
    pkhbt  r12, r4, lr, lsl #16     // T2
    // ---- T3 = mont(OMEGA * (T1 - T2)) ----
    ssub16 r4, r11, r12             // D = T1 - T2
    movw   r0, #64650               // OMEGA (low16)
    smulbt lr, r0, r4               // OMEGA*D.hi
    smulbb r4, r0, r4               // OMEGA*D.lo
    mul    r0, r4, r2
    sxth   r0, r0
    smulbb r0, r0, r3
    sub    r4, r4, r0
    asr    r4, r4, #16              // t3.lo
    mul    r0, lr, r2               // (r0 was OMEGA, freed after smulbb above)
    sxth   r0, r0
    smulbb r0, r0, r3
    sub    lr, lr, r0
    asr    lr, lr, #16              // t3.hi
    pkhbt  r0, r4, lr, lsl #16      // T3  (in r0)
    // ---- combines (A loaded into r4) ----
    ldr    r4, [r10]                // A = r[i]|r[i+1]
    ssub16 lr, r4, r11              // A - T1
    ssub16 lr, lr, r0               //   - T3   -> r[i+2step]
    str    lr, [r9], #4
    ssub16 lr, r4, r12              // A - T2
    sadd16 lr, lr, r0               //   + T3   -> r[i+step]
    str    lr, [r8], #4
    sadd16 lr, r4, r11              // A + T1
    sadd16 lr, lr, r12              //   + T2   -> r[i]
    str    lr, [r10], #4
    subs   r5, r5, #1
    bne    .Lr3_bfly
    // next start: start_byte += 3*step_bytes
    ldr    r4, [sp, #8]
    ldr    r11, [sp, #16]
    add    r11, r11, r4
    str    r11, [sp, #16]
    ldr    r12, [sp, #12]           // n_bytes
    cmp    r11, r12
    blt    .Lr3_start
    add    sp, sp, #20
    pop    {r4-r11, pc}
.size ntt_radix3_layer, .-ntt_radix3_layer
