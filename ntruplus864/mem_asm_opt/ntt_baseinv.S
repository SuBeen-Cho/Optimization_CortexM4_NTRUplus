// SPDX-License-Identifier: MIT
// NTRU+864 (X^3-zeta, degree-3) Cortex-M4 baseinv: keygen-only base-ring inversion.
// Bit-identical to the C baseinv + fqinv (host-verified in results/ntruplus_opt/
// baseinv_check.c, 5M random incl. the non-invertible early-out); KAT-gated.
// Register-resident a0..a2 / r0..r2, inlined fqinv addition chain, short Montgomery.
//   int baseinv(int16_t r[3], const int16_t a[3], int16_t zeta);
// Returns 0 (invertible, r[0..2] written) or 1 (non-invertible, r[] untouched).
// QINV=12929 Q=3457 RINV=-682 (u16 64854)

.syntax unified
.cpu cortex-m4
.thumb

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
    ldrsh  r7, [r1, #4]             // a2   (r1 free; r2=zeta)
    movw   r3, #12929               // QINV
    movw   r4, #3457                // Q
    // r0v = mont(a1*a2)
    smulbb lr, r6, r7
    BMONT  r8, lr                   // r0v
    // r1v = mont(a2*a2)
    smulbb lr, r7, r7
    BMONT  r9, lr                   // r1v
    // r2v = mont(a1*a1 - a0*a2)
    smulbb lr, r6, r6
    smulbb r1, r5, r7
    sub    lr, lr, r1
    BMONT  r10, lr                  // r2v
    // r0v = mont(a0*a0 - r0v*zeta)
    smulbb lr, r5, r5
    smulbb r1, r8, r2
    sub    lr, lr, r1
    BMONT  r8, lr                   // r0v
    // r1v = mont(r1v*zeta - a0*a1)
    smulbb lr, r9, r2
    smulbb r1, r5, r6
    sub    lr, lr, r1
    BMONT  r9, lr                   // r1v
    // t = mont(r2v*a1 + r1v*a2)
    smulbb lr, r10, r6
    smlabb lr, r9, r7, lr
    BMONT  r11, lr                  // t
    // t = mont(t*zeta + r0v*a0)
    smulbb lr, r11, r2
    smlabb lr, r8, r5, lr
    BMONT  r11, lr                  // t
    cmp    r11, #0
    beq    .Lnoninv
    // t = fqinv(t) ; a=r11, locals t1=r5 t2=r6 t3l=r7 (a0..a2 dead), tinv->r6
    smulbb lr, r11, r11
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
    smulbb lr, r6, r11
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
    // [FIX-ASM] clean convention: clean's fqinv returns t2 directly. Removed the
    // trailing mont(RINV*t2) that divided by R (asm baseinv = clean*R^-1). r6 now
    // holds clean's fqinv(t) so final scaling matches clean.
    // final scaling (no negation for deg-3)
    smulbb lr, r8, r6
    BMONT  r5, lr
    strh   r5, [r0]                 // r[0]
    smulbb lr, r9, r6
    BMONT  r5, lr
    strh   r5, [r0, #2]             // r[1]
    smulbb lr, r10, r6
    BMONT  r5, lr
    strh   r5, [r0, #4]             // r[2]
    mov    r0, #0
    pop    {r4-r11, pc}
.Lnoninv:
    mov    r0, #1
    pop    {r4-r11, pc}
.size baseinv, .-baseinv
