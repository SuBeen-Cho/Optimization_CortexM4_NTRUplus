// SPDX-License-Identifier: MIT
// NTRU+768 Cortex-M4 NTT assembly (q=3457, Montgomery R=2^16, QINV=12929).
//
// ntt_ct_layers(): radix-2 Cooley-Tukey butterfly layers (step 64..4) + the final
// Barrett pass of ntt(). PACKED-16: within each block the twiddle zeta is constant,
// so two adjacent butterflies (i,i+1) share it. We load r[i]|r[i+1] and
// r[i+step]|r[i+step+1] as 32-bit words, Montgomery-multiply both halves of the
// second word by zeta (SMULBB low / SMULBT high), repack with PKHBT, then do the
// add/sub of BOTH butterflies in one SADD16 / SSUB16. Halves load/store + add/sub
// traffic vs the scalar version. Bit-identical to the C reference; KAT-gated.
//
//   void ntt_ct_layers(int16_t r[768], const int16_t *zp);  // zp = &zetas[6]
//
// montgomery_reduce(p): m=(int16)(p*QINV); return (p - m*q) >> 16;  (no reduce on +/-)
// barrett_reduce(a):    t=(v*a + (1<<25))>>26; return a - t*q;  v=19412

.syntax unified
.cpu cortex-m4
.thumb

.global ntt_ct_layers
.type   ntt_ct_layers, %function
.align 2
ntt_ct_layers:
    push   {r4-r11, lr}
    sub    sp, sp, #4
    str    r0, [sp]                 // base r
    movw   r2, #12929               // qinv
    movw   r3, #3457                // q
    mov    r4, #64                  // step (coeffs)
.Llayer:
    ldr    r8, [sp]                 // pA = base
    lsl    r5, r4, #1               // step_bytes = step*2
    movw   r6, #384
    udiv   r6, r6, r4               // nblocks = 384/step
.Lblock:
    ldrsh  r7, [r1], #2             // zeta = *zp++
    add    r9, r8, r5               // pB = pA + step_bytes
    asr    r10, r4, #1              // cnt = step/2  (packed butterflies per block)
.Lbfly:
    ldr    r11, [r8]                // A = r[i] | r[i+1]
    ldr    r12, [r9]                // B = r[i+step] | r[i+step+1]
    smulbt lr,  r7, r12             // hi prod = zeta * r[i+step+1]
    smulbb r12, r7, r12             // lo prod = zeta * r[i+step]
    // montgomery_reduce(lo prod) -> r12
    mul    r0,  r12, r2
    sxth   r0,  r0
    smulbb r0,  r0,  r3
    sub    r12, r12, r0
    asr    r12, r12, #16            // t_lo
    // montgomery_reduce(hi prod) -> lr
    mul    r0,  lr,  r2
    sxth   r0,  r0
    smulbb r0,  r0,  r3
    sub    lr,  lr,  r0
    asr    lr,  lr,  #16            // t_hi
    pkhbt  r12, r12, lr, lsl #16    // T = t_lo | t_hi<<16
    sadd16 lr,  r11, r12            // r[i]|r[i+1]        = A + T
    ssub16 r11, r11, r12            // r[i+step]|r[i+..+1]= A - T
    str    lr,  [r8], #4
    str    r11, [r9], #4
    subs   r10, r10, #1
    bne    .Lbfly
    mov    r8, r9                   // pA = pB (next block start)
    subs   r6, r6, #1
    bne    .Lblock
    asr    r4, r4, #1               // step >>= 1
    cmp    r4, #4
    bge    .Llayer

    // Barrett pass over all 768 coeffs, packed 2-at-a-time:
    //   a -= ((v*a + (1<<25))>>26)*q
    ldr    r8, [sp]                 // p = base
    movw   r4, #19412               // v
    add    r10, r8, #(768*2)        // end
.Lbar:
    ldr    r6,  [r8]                // A = a_lo | a_hi
    sxth   r11, r6                  // a_lo (signed)
    asr    r12, r6, #16             // a_hi (signed)
    smulbb r9,  r4, r6              // v*a_lo
    add    r9,  r9, #(1<<25)
    asr    r9,  r9, #26
    smulbb r9,  r9, r3              // t_lo*q
    sub    r11, r11, r9             // a_lo - t_lo*q
    smulbt r9,  r4, r6              // v*a_hi
    add    r9,  r9, #(1<<25)
    asr    r9,  r9, #26
    smulbb r9,  r9, r3              // t_hi*q
    sub    r12, r12, r9             // a_hi - t_hi*q
    pkhbt  r11, r11, r12, lsl #16
    str    r11, [r8], #4
    cmp    r8, r10
    bne    .Lbar

    add    sp, sp, #4
    pop    {r4-r11, pc}
.size ntt_ct_layers, .-ntt_ct_layers
