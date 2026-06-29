// SPDX-License-Identifier: MIT

#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

typedef struct{
	int16_t coeffs[NTRUPLUS_N];
} poly;

void poly_tobytes(uint8_t r[NTRUPLUS_POLYBYTES], const poly *a);
void poly_frombytes(poly *r, const uint8_t a[NTRUPLUS_POLYBYTES]);

void poly_cbd1(poly *r, const uint8_t buf[NTRUPLUS_N/4]);
void poly_sotp(poly *r, const unsigned char *msg, const unsigned char *buf);
int  poly_sotp_inv(unsigned char *msg, const poly *e, const unsigned char *buf);

void poly_ntt(poly *r, const poly *a);
void poly_ntt_ternary(poly *r, const poly *a);
void poly_invntt(poly *r, const poly *a);
int  poly_baseinv(poly *r, const poly *a);
void poly_basemul(poly *r, const poly *a, const poly *b);
void poly_basemul_add(poly *r, const poly *a, const poly *b, const poly *c);
void poly_sub(poly *r, const poly *a, const poly *b);
void poly_triple(poly *r, const poly *a);
void poly_crepmod3(poly *r, const poly *a);

/* [OPT6] Fused INVNTT(lazy Barrett) + crepmod3 */
void poly_invntt_crepmod3(poly *r, const poly *a);

/* [OPT7] Ternary compression */
void poly_compress_ternary(uint8_t r[NTRUPLUS_N/4], const poly *a);
void poly_decompress_ternary(poly *r, const uint8_t a[NTRUPLUS_N/4]);

/* [OPT9] Streaming CBD1 from seed */
void poly_cbd1_from_seed(poly *r, const uint8_t seed[32]);

/* [OPT7] basemul with one operand from tobytes bytes (eliminates poly b in keygen) */
void poly_basemul_packed_tobytes(poly *r, const poly *a, const uint8_t *b_tobytes);

/* [OPT8] Sequential 12-bit packing */
void poly_pack_seq(uint8_t r[NTRUPLUS_POLYBYTES], const poly *a);
void poly_basemul_add_packed(poly *r, const poly *a, const uint8_t *b_packed, const poly *c);

/* [OPT5] Both operands from packed bytes (eliminates poly b in encaps) */
void poly_basemul_add_2packed(poly *r, const uint8_t *a_tobytes,
                              const uint8_t *b_packseq, const poly *c);

/* [OPT3] Fused serialize+hash_g */
void poly_tobytes_hash_g(uint8_t *hash_out, const poly *a);

/* [OPT3] Constant-time polynomial comparison */
int poly_verify(const poly *a, const poly *b);

/* [OPT8] Decaps streaming functions */
void poly_basemul_2tobytes(poly *r, const uint8_t *a_tobytes, const uint8_t *b_tobytes);
void poly_sub_basemul_tobytes(poly *r, const uint8_t *ct_tobytes,
                               const poly *a, const uint8_t *hinv_tobytes);
int poly_sotp_inv_compressed(unsigned char *msg, const uint8_t *m1t, const unsigned char *buf);

#endif
