#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

extern const int16_t zetas[288];

void ntt(int16_t r[NTRUPLUS_N], const int16_t a[NTRUPLUS_N]);
void invntt(int16_t r[NTRUPLUS_N], const int16_t a[NTRUPLUS_N]);

// asm: radix-2 CT layers (barrett-inside) of ntt(); zp = &zetas[18]
void ntt_ct_layers(int16_t *r, const int16_t *zp);
// asm: one radix-3 (OMEGA) layer of ntt()
void ntt_radix3_layer(int16_t *r, const int16_t *zp, int n, int step);
void invntt_gs_layer(int16_t *r, const int16_t *zp, int n, int step);
void invntt_scale(int16_t *r, int n, int ninv, int twoninv);
void invntt_radix3_layer(int16_t *r, const int16_t *zp, int n, int step, int do_barrett);

int  baseinv(int16_t r[3], const int16_t a[3], int16_t zeta);
void basemul(int16_t r[3], const int16_t a[3], const int16_t b[3], int16_t zeta);
void basemul_add(int16_t r[3], const int16_t a[3], const int16_t b[3], const int16_t c[3], int16_t zeta);

#endif
