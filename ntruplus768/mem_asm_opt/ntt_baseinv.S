// SPDX-License-Identifier: MIT
// NTRU+ (X^4-zeta, degree-4) Cortex-M4 baseinv: keygen-only base-ring inversion.
// Bit-identical to the C baseinv + fqinv (host-verified in results/ntruplus_opt/
// baseinv_check.c, 5M random incl. the non-invertible early-out); KAT-gated.
// Register-resident a0..a3 / t0..t3, inlined fqinv addition chain, short Montgomery.
//   int baseinv(int16_t r[4], const int16_t a[4], int16_t zeta);
// Returns 0 (invertible, r[0..3] written) or 1 (non-invertible, r[] untouched —
// poly_baseinv zeroes the poly on 1).   QINV=12929 Q=3457 RINV=-682 (u16 64854)

.syntax unified
.cpu cortex-m4
.thumb

// montgomery_reduce of the 32-bit accumulator \acc -> \dst (int16).
// m=(int16)(acc*QINV); dst=(acc - m*Q)>>16.  r3=QINV, r4=Q, r1=scratch.
.macro BMONT dst, acc
    smulbb r1, \acc, r3
    smulbb r1, r1, r4
    sub    \acc, \acc, r1
    asr    \dst, \acc, #16
.endm

.global baseinv
.type   baseinv, %function
.align 2
baseinv:
    push   {r4-r11, lr}
    ldrsh  r5, [r1]                 // a0
    ldrsh  r6, [r1, #2]             // a1
    ldrsh  r7, [r1, #4]             // a2
    ldrsh  r8, [r1, #6]             // a3   (r1 now free as BMONT scratch; r2=zeta)
    movw   r3, #12929               // QINV
    movw   r4, #3457                // Q
    // t0 = mont(a2*a2 - 2*a1*a3)
    smulbb lr, r7, r7
    smulbb r1, r6, r8
    sub    lr, lr, r1, lsl #1
    BMONT  r9, lr                   // t0
    // t1 = mont(a3*a3)
    smulbb lr, r8, r8
    BMONT  r10, lr                  // t1
    // t0 = mont(a0*a0 + t0*zeta)
    smulbb lr, r5, r5
    smlabb lr, r9, r2, lr
    BMONT  r9, lr                   // t0
    // t1 = mont(a1*a1 + t1*zeta - 2*a0*a2)
    smulbb lr, r6, r6
    smlabb lr, r10, r2, lr
    smulbb r1, r5, r7
    sub    lr, lr, r1, lsl #1
    BMONT  r10, lr                  // t1
    // t2 = mont(t1*zeta)
    smulbb lr, r10, r2
    BMONT  r11, lr                  // t2
    // t3 = mont(t0*t0 - t1*t2)
    smulbb lr, r9, r9
    smulbb r1, r10, r11
    sub    lr, lr, r1
    BMONT  r12, lr                  // t3
    cmp    r12, #0
    beq    .Lnoninv
    // r[0..3] = mont(...) ; store immediately (frees a0..a3 / t0..t2 after)
    smulbb lr, r5, r9               // a0*t0
    smlabb lr, r7, r11, lr          // + a2*t2
    BMONT  lr, lr
    strh   lr, [r0]                 // r[0]
    smulbb lr, r8, r11              // a3*t2
    smlabb lr, r6, r9, lr           // + a1*t0
    BMONT  lr, lr
    strh   lr, [r0, #2]             // r[1]
    smulbb lr, r7, r9               // a2*t0
    smlabb lr, r5, r10, lr          // + a0*t1
    BMONT  lr, lr
    strh   lr, [r0, #4]             // r[2]
    smulbb lr, r6, r10              // a1*t1
    smlabb lr, r8, r9, lr           // + a3*t0
    BMONT  lr, lr
    strh   lr, [r0, #6]             // r[3]
    // t3 = fqinv(t3) ; a=r12, locals t1=r5 t2=r6 t3l=r7 (a0..a3 dead), tinv->r6
    smulbb lr, r12, r12
    BMONT  r5, lr                   // t1 = a*a
    smulbb lr, r5, r5
    BMONT  r6, lr                   // t2 = t1*t1
    smulbb lr, r6, r6
    BMONT  r6, lr                   // t2 = t2*t2
    smulbb lr, r6, r6
    BMONT  r7, lr                   // t3 = t2*t2
    smulbb lr, r5, r6
    BMONT  r5, lr                   // t1 = t1*t2
    smulbb lr, r5, r7
    BMONT  r6, lr                   // t2 = t1*t3
    smulbb lr, r6, r6
    BMONT  r6, lr                   // t2 = t2*t2
    smulbb lr, r6, r12
    BMONT  r6, lr                   // t2 = t2*a
    smulbb lr, r5, r6
    BMONT  r5, lr                   // t1 = t1*t2
    smulbb lr, r6, r6
    BMONT  r6, lr                   // t2 = t2*t2  (x6)
    smulbb lr, r6, r6
    BMONT  r6, lr
    smulbb lr, r6, r6
    BMONT  r6, lr
    smulbb lr, r6, r6
    BMONT  r6, lr
    smulbb lr, r6, r6
    BMONT  r6, lr
    smulbb lr, r6, r6
    BMONT  r6, lr
    smulbb lr, r6, r5
    BMONT  r6, lr                   // t2 = t2*t1  == clean fqinv output (tinv)
    // [FIX-ASM] clean convention: clean's fqinv returns t2 directly (no extra
    // Montgomery factor). The previous trailing mont(RINV*t2) divided by R,
    // making asm baseinv = clean * R^-1 and diverging from clean serialization.
    // Removed — r6 now holds clean's fqinv(t3), so final scaling matches clean.
    // final scaling: r[i] = (+/-) mont(r[i]*tinv)
    ldrsh  r5, [r0]
    smulbb lr, r5, r6
    BMONT  lr, lr
    strh   lr, [r0]                 // r[0]
    ldrsh  r5, [r0, #2]
    smulbb lr, r5, r6
    BMONT  lr, lr
    rsb    lr, lr, #0               // r[1] = -mont(...)
    strh   lr, [r0, #2]
    ldrsh  r5, [r0, #4]
    smulbb lr, r5, r6
    BMONT  lr, lr
    strh   lr, [r0, #4]             // r[2]
    ldrsh  r5, [r0, #6]
    smulbb lr, r5, r6
    BMONT  lr, lr
    rsb    lr, lr, #0               // r[3] = -mont(...)
    strh   lr, [r0, #6]
    mov    r0, #0
    pop    {r4-r11, pc}
.Lnoninv:
    mov    r0, #1
    pop    {r4-r11, pc}
.size baseinv, .-baseinv
