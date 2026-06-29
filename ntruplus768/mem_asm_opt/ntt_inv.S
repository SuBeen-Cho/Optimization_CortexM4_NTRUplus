// SPDX-License-Identifier: MIT
// NTRU+ Cortex-M4 invntt() helpers (q=3457). Bit-identical to the C; KAT-gated.
// C drives the step-loops and the (decreasing) zeta index; asm does each layer.
//   invntt_gs_layer(r, zp, n, step): one Gentleman-Sande radix-2 layer.
//     zeta=*zp-- per start; per i: t1=r[i+step];
//       r[i+step]=mont(zeta, t1-r[i]);  r[i]=barrett(r[i]+t1)
//   invntt_scale(r, n, ninv, twoninv): final X^(n/2) merge + scaling.
//     t1=r[i]+r[i+n/2]; t2=mont(ZM5, r[i]-r[i+n/2]);
//     r[i]=mont(ninv, t1-t2); r[i+n/2]=mont(twoninv, t2)
// QINV=12929 Q=3457 ZM5=-1665(=63871) BARV=19412

.syntax unified
.cpu cortex-m4
.thumb

// ---------------- invntt_gs_layer(r0=r, r1=zp, r2=n, r3=step) ----------------
.global invntt_gs_layer
.type   invntt_gs_layer, %function
.align 2
// PACKED-16: two adjacent GS butterflies (i,i+1) share zeta within a block.
//   D = B - A (ssub16);  r[i+s]|r[i+s+1] = mont(zeta*D.lo) | mont(zeta*D.hi)
//   S = A + B (sadd16);  r[i]|r[i+1]     = barrett(S.lo)   | barrett(S.hi)
// Bit-identical to the scalar GS (all ops are per-lane int16, same as the C).
invntt_gs_layer:
    push   {r4-r11, lr}
    sub    sp, sp, #12
    lsl    r2, r2, #1               // n_bytes
    str    r2, [sp, #0]
    lsl    r3, r3, #1               // step_bytes
    str    r3, [sp, #4]
    mov    r6, #0
    str    r6, [sp, #8]             // start_byte
    movw   r2, #12929               // qinv
    movw   r3, #3457                // q
    movw   r4, #19412               // v (BARV)
.Lgs_start:
    ldr    r6, [sp, #8]             // start_byte
    add    r6, r0, r6               // pA = base + start_byte
    ldr    r9, [sp, #4]             // step_bytes
    add    r7, r6, r9               // pB = pA + step_bytes
    asr    r8, r9, #2               // cnt = step_bytes/4 = step/2
    ldrsh  r5, [r1], #-2            // zeta = *zp--
.Lgs_bfly:
    ldr    r9,  [r6]                // A = r[i] | r[i+1]
    ldr    r10, [r7]                // B = r[i+step] | r[i+step+1]
    ssub16 r11, r10, r9             // D = B - A  (packed diffs)
    smulbt r12, r5, r11             // zeta * D.hi
    smulbb r11, r5, r11             // zeta * D.lo
    mul    lr,  r11, r2             // mont(D.lo)
    sxth   lr,  lr
    smulbb lr,  lr,  r3
    sub    r11, r11, lr
    asr    r11, r11, #16            // t_lo
    mul    lr,  r12, r2             // mont(D.hi)
    sxth   lr,  lr
    smulbb lr,  lr,  r3
    sub    r12, r12, lr
    asr    r12, r12, #16            // t_hi
    pkhbt  r11, r11, r12, lsl #16
    str    r11, [r7], #4            // r[i+step]|r[i+step+1]
    sadd16 r11, r9, r10             // S = A + B  (packed sums)
    sxth   r9,  r11                 // s_lo
    asr    r10, r11, #16            // s_hi
    smulbb r12, r4, r11             // v*s_lo
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3             // t*q
    sub    r9,  r9,  r12            // barrett(s_lo)
    smulbt r12, r4, r11             // v*s_hi
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    sub    r10, r10, r12            // barrett(s_hi)
    pkhbt  r9,  r9,  r10, lsl #16
    str    r9,  [r6], #4            // r[i]|r[i+1]
    subs   r8,  r8,  #1
    bne    .Lgs_bfly
    // odd-step tail (step=3 for N=864): one lone scalar butterfly (pA,pB now aligned)
    ldr    r9, [sp, #4]             // step_bytes
    tst    r9, #2                   // step odd?  (step_bytes=2*step; bit1 set iff step odd)
    beq    .Lgs_notail
    ldrsh  r10, [r7]                // t1 = r[i+step]
    ldrsh  r11, [r6]                // a = r[i]
    sub    r12, r10, r11            // t1 - a
    smulbb r12, r5, r12             // zeta*(t1-a)
    mul    lr, r12, r2
    sxth   lr, lr
    smulbb lr, lr, r3
    sub    r12, r12, lr
    asr    r12, r12, #16            // mont -> r[i+step]
    strh   r12, [r7]
    add    r12, r11, r10            // a + t1
    smulbb lr, r4, r12              // v*(a+t1)
    add    lr, lr, #(1<<25)
    asr    lr, lr, #26
    smulbb lr, lr, r3
    sub    r12, r12, lr             // barrett -> r[i]
    strh   r12, [r6]
.Lgs_notail:
    ldr    r6, [sp, #8]
    ldr    r9, [sp, #4]
    add    r6, r6, r9, lsl #1       // start += 2*step_bytes
    str    r6, [sp, #8]
    ldr    r10, [sp, #0]
    cmp    r6, r10
    blt    .Lgs_start
    add    sp, sp, #12
    pop    {r4-r11, pc}
.size invntt_gs_layer, .-invntt_gs_layer

// -------------- invntt_scale(r0=r, r1=n, r2=ninv, r3=twoninv) --------------
.global invntt_scale
.type   invntt_scale, %function
.align 2
// PACKED-16: two adjacent pairs (i,i+1) of the X^(n/2) merge at once.
//   Dm=A-B (ssub16); T2=mont(ZM5*Dm) packed; T1=A+B (sadd16);
//   r[i]=mont(ninv*(T1-T2)); r[i+n/2]=mont(twoninv*T2)   -- all per-lane, == C.
invntt_scale:
    push   {r4-r11, lr}
    add    r4, r0, r1               // pB = &r[n/2]  (n bytes)
    lsr    r5, r1, #2               // cnt = (n/2)/2 = n/4
    movw   r6, #12929               // qinv
    movw   r7, #3457                // q
    movw   r8, #63871               // ZM5 = -1665
.Lsc:
    ldr    r9,  [r0]                // A = r[i]|r[i+1]
    ldr    r10, [r4]                // B = r[i+n/2]|..
    ssub16 r11, r9, r10             // Dm = A - B
    smulbt r12, r8, r11             // ZM5 * Dm.hi
    smulbb r11, r8, r11             // ZM5 * Dm.lo
    mul    lr, r11, r6              // mont(Dm.lo)
    sxth   lr, lr
    smulbb lr, lr, r7
    sub    r11, r11, lr
    asr    r11, r11, #16            // t2_lo
    mul    lr, r12, r6              // mont(Dm.hi)
    sxth   lr, lr
    smulbb lr, lr, r7
    sub    r12, r12, lr
    asr    r12, r12, #16            // t2_hi
    pkhbt  r11, r11, r12, lsl #16   // T2 = t2_lo|t2_hi
    sadd16 r12, r9, r10             // T1 = A + B
    ssub16 r9,  r12, r11            // T1 - T2  (packed)
    smulbt r10, r2, r9              // ninv * (T1-T2).hi
    smulbb r9,  r2, r9              // ninv * (T1-T2).lo
    mul    lr, r9, r6               // mont lo
    sxth   lr, lr
    smulbb lr, lr, r7
    sub    r9, r9, lr
    asr    r9, r9, #16
    mul    lr, r10, r6              // mont hi
    sxth   lr, lr
    smulbb lr, lr, r7
    sub    r10, r10, lr
    asr    r10, r10, #16
    pkhbt  r9, r9, r10, lsl #16
    str    r9, [r0], #4             // r[i]|r[i+1]
    smulbt r10, r3, r11             // twoninv * T2.hi
    smulbb r11, r3, r11             // twoninv * T2.lo
    mul    lr, r11, r6              // mont lo
    sxth   lr, lr
    smulbb lr, lr, r7
    sub    r11, r11, lr
    asr    r11, r11, #16
    mul    lr, r10, r6              // mont hi
    sxth   lr, lr
    smulbb lr, lr, r7
    sub    r10, r10, lr
    asr    r10, r10, #16
    pkhbt  r11, r11, r10, lsl #16
    str    r11, [r4], #4            // r[i+n/2]|..
    subs   r5, r5, #1
    bne    .Lsc
    pop    {r4-r11, pc}
.size invntt_scale, .-invntt_scale

// -- invntt_radix3_layer(r0=r, r1=zp, r2=n, r3=step, [sp+36]=do_barrett) --
// zeta2=zp[0], zeta1=zp[-1], zp-=2 per start; per i:
//   t1=mont(OMEGA,b-a); t2=mont(zeta1,c-a+t1); t3=mont(zeta2,c-b-t1);
//   r[i]=(barrett?)(a+b+c); r[i+step]=t2; r[i+2step]=t3   (a=r[i],b=r[i+s],c=r[i+2s])
.global invntt_radix3_layer
.type   invntt_radix3_layer, %function
.align 2
invntt_radix3_layer:
    push   {r4-r11, lr}
    ldr    r4, [sp, #36]            // do_barrett
    sub    sp, sp, #28
    str    r0, [sp, #24]            // rbase
    lsl    r2, r2, #1               // n_bytes
    str    r2, [sp, #0]
    lsl    r3, r3, #1               // sb = step_bytes
    str    r3, [sp, #4]
    add    r5, r3, r3, lsl #1       // 3*sb
    str    r5, [sp, #8]
    mov    r5, #0
    str    r5, [sp, #12]            // start_byte
    str    r4, [sp, #16]            // do_barrett
    str    r1, [sp, #20]            // zp
    movw   r2, #12929               // qinv
    movw   r3, #3457                // q
.Lir3_start:
    ldr    r1, [sp, #20]            // zp
    ldrsh  r7, [r1]                 // zeta2 = zp[0]
    ldrsh  r6, [r1, #-2]            // zeta1 = zp[-1]
    sub    r1, r1, #4
    str    r1, [sp, #20]            // zp -= 2 elements
    ldr    r0, [sp, #24]            // rbase
    ldr    r1, [sp, #12]            // start_byte
    add    r0, r0, r1               // pA = &r[start]
    ldr    r4, [sp, #4]             // sb (step_bytes); pC offset = [r0, r4, lsl#1]
    add    r1, r0, r4               // pA_end = pA + sb
.Lir3_bfly:
    // PACKED-16: two adjacent triples (i,i+1) share (zeta1,zeta2,OMEGA).
    ldr    r8,  [r0]                // A = r[i]|r[i+1]
    ldr    r9,  [r0, r4]            // B = r[i+s]|r[i+s+1]
    ldr    r10, [r0, r4, lsl #1]    // C = r[i+2s]|r[i+2s+1]
    // T1 = mont(OMEGA * (B - A))
    ssub16 r12, r9, r8              // B - A
    movw   r5, #64650               // OMEGA (low16)
    smulbt lr, r5, r12              // OMEGA*hi
    smulbb r12, r5, r12             // OMEGA*lo
    mul    r5, r12, r2
    sxth   r5, r5
    smulbb r5, r5, r3
    sub    r12, r12, r5
    asr    r12, r12, #16            // t1.lo
    mul    r5, lr, r2
    sxth   r5, r5
    smulbb r5, r5, r3
    sub    lr, lr, r5
    asr    lr, lr, #16              // t1.hi
    pkhbt  r11, r12, lr, lsl #16    // T1
    // T2 = mont(zeta1 * ((C - A) + T1))  -> r[i+s]
    ssub16 r12, r10, r8             // C - A
    sadd16 r12, r12, r11            //   + T1
    smulbt lr, r6, r12              // z1*hi
    smulbb r12, r6, r12             // z1*lo
    mul    r5, r12, r2
    sxth   r5, r5
    smulbb r5, r5, r3
    sub    r12, r12, r5
    asr    r12, r12, #16            // t2.lo
    mul    r5, lr, r2
    sxth   r5, r5
    smulbb r5, r5, r3
    sub    lr, lr, r5
    asr    lr, lr, #16              // t2.hi
    pkhbt  r12, r12, lr, lsl #16    // T2
    str    r12, [r0, r4]            // r[i+s] = T2
    // T3 = mont(zeta2 * ((C - B) - T1))  -> r[i+2s]
    ssub16 r12, r10, r9             // C - B
    ssub16 r12, r12, r11            //   - T1
    smulbt lr, r7, r12              // z2*hi
    smulbb r12, r7, r12             // z2*lo
    mul    r5, r12, r2
    sxth   r5, r5
    smulbb r5, r5, r3
    sub    r12, r12, r5
    asr    r12, r12, #16            // t3.lo
    mul    r5, lr, r2
    sxth   r5, r5
    smulbb r5, r5, r3
    sub    lr, lr, r5
    asr    lr, lr, #16              // t3.hi
    pkhbt  r12, r12, lr, lsl #16    // T3
    str    r12, [r0, r4, lsl #1]    // r[i+2s] = T3
    // S = A + B + C ; optional per-lane barrett -> r[i]
    sadd16 r8, r8, r9               // A + B
    sadd16 r8, r8, r10              //   + C  (S; B,C now dead)
    ldr    r5, [sp, #16]            // do_barrett
    cmp    r5, #0
    beq    .Lir3_nob
    movw   r5, #19412               // BARV v
    sxth   r9, r8                   // S.lo
    smulbb r12, r5, r9
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    sub    r9, r9, r12              // barrett(S.lo)
    asr    r10, r8, #16             // S.hi
    smulbb r12, r5, r10
    add    r12, r12, #(1<<25)
    asr    r12, r12, #26
    smulbb r12, r12, r3
    sub    r10, r10, r12            // barrett(S.hi)
    pkhbt  r8, r9, r10, lsl #16
.Lir3_nob:
    str    r8, [r0], #4             // r[i]|r[i+1], pA += 2 coeffs
    cmp    r0, r1
    blt    .Lir3_bfly
    ldr    r4, [sp, #12]            // start_byte
    ldr    r5, [sp, #8]             // 3*sb
    add    r4, r4, r5
    str    r4, [sp, #12]            // start += 3*sb
    ldr    r6, [sp, #0]             // n_bytes
    cmp    r4, r6
    blt    .Lir3_start
    add    sp, sp, #28
    pop    {r4-r11, pc}
.size invntt_radix3_layer, .-invntt_radix3_layer
