// SPDX-License-Identifier: MIT

#include <stddef.h>
#include <stdint.h>
#include "api.h"
#include "params.h"
#include "symmetric.h"
#include "poly.h"
#include "verify.h"
#include "fips202.h"
#include "randombytes.h"

/*************************************************
* [OPT1] poly 5->3, [OPT4] poly 3->2, [OPT9] streaming CBD1
* [OPT7] poly 2->1: in-place baseinv + basemul from sk bytes
**************************************************/
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{
	uint8_t seed[32];

	poly a;  /* [OPT7] Only 1 poly! b removed via in-place baseinv */

	/* Phase 1: f and finv */
	do {
		randombytes(seed, 32);
		poly_cbd1_from_seed(&a, seed);
		poly_triple(&a, &a);
		a.coeffs[0] += 1;
		poly_ntt(&a, &a);
		poly_tobytes(sk, &a);          /* [OPT7] Save f BEFORE baseinv */
	} while(poly_baseinv(&a, &a));     /* [OPT7] In-place baseinv! */

	poly_tobytes(sk + NTRUPLUS_POLYBYTES, &a);  /* Save finv */

	/* Phase 2: g, h, hinv */
	do {
		randombytes(seed, 32);
		poly_cbd1_from_seed(&a, seed);
		poly_triple(&a, &a);
		poly_ntt(&a, &a);

		/* [OPT7] Read finv directly from sk bytes, no poly b needed */
		poly_basemul_packed_tobytes(&a, &a, sk + NTRUPLUS_POLYBYTES);
		poly_tobytes(pk, &a);          /* [OPT7] Save h BEFORE baseinv */
	} while(poly_baseinv(&a, &a));     /* [OPT7] In-place baseinv! */

	poly_tobytes(sk + NTRUPLUS_POLYBYTES, &a);  /* Save hinv */
	hash_f(sk + 2 * NTRUPLUS_POLYBYTES, pk);

	return 0;
}

/*************************************************
* [OPT1] poly 4->3, [OPT4] buf2->hashbuf, [OPT8] poly 3->2
* [OPT5] poly 2->1 (eliminate b via 2packed), [OPT10] union,
* [OPT11] direct ss output
**************************************************/
int crypto_kem_enc(unsigned char *ct,
                   unsigned char *ss,
                   const unsigned char *pk)
{
	uint8_t msg[NTRUPLUS_N / 8 + NTRUPLUS_SYMBYTES];

	/* [OPT10] Union: cbd_seed and hashbuf share memory */
	union {
		uint8_t cbd_seed[NTRUPLUS_N / 4];
		uint8_t hashbuf[NTRUPLUS_N / 4];
	} enc_u;

	poly a;  /* [OPT5] poly b removed: read pk directly from bytes */

	randombytes(msg, NTRUPLUS_N / 8);
	hash_f(msg + NTRUPLUS_N / 8, pk);

	/* [OPT11] Squeeze ss directly to output, cbd_seed to union */
	{
		shake256incctx hctx;
		uint8_t domain = 0x2;
		shake256_inc_init(&hctx);
		shake256_inc_absorb(&hctx, &domain, 1);
		shake256_inc_absorb(&hctx, msg, NTRUPLUS_N / 8 + NTRUPLUS_SYMBYTES);
		shake256_inc_finalize(&hctx);
		shake256_inc_squeeze(ss, NTRUPLUS_SSBYTES, &hctx);
		shake256_inc_squeeze(enc_u.cbd_seed, NTRUPLUS_N / 4, &hctx);
		shake256_inc_ctx_release(&hctx);
	}

	poly_cbd1(&a, enc_u.cbd_seed);
	poly_ntt_ternary(&a, &a);  /* [OPT-NTT] cbd1 output is ternary */

	/* [OPT3] Fused serialize+hash */
	poly_tobytes_hash_g(enc_u.hashbuf, &a);

	/* [OPT8] Pack r into ct scratch, then streaming basemul_add */
	poly_pack_seq(ct, &a);

	poly_sotp(&a, msg, enc_u.hashbuf);
	poly_ntt(&a, &a);

	/* [OPT5] Read h from pk (tobytes) and r from ct (pack_seq) directly */
	poly_basemul_add_2packed(&a, pk, ct, &a);

	poly_tobytes(ct, &a);

	return 0;
}

/*************************************************
* [OPT8] Decaps 1-poly: hash-based verify + streaming basemul
* poly 2->1 via hash_g comparison instead of poly_verify
**************************************************/
int crypto_kem_dec(unsigned char *ss,
                   const unsigned char *ct,
                   const unsigned char *sk)
{
	uint8_t msg[NTRUPLUS_N / 8 + NTRUPLUS_SYMBYTES];

	/* [OPT8] r2_hash stored separately for verify */
	uint8_t r2_hash[NTRUPLUS_N / 4];

	/* [OPT9] m1_ternary folded into union with buf3: after dedup, hashbuf is
	 * gone (sotp_inv now uses r2_hash), so m1_ternary (live steps 4-7) and
	 * buf3 (live step 9+) have disjoint lifetimes and can share storage. */
	union {
		uint8_t m1_ternary[NTRUPLUS_N / 4];
		uint8_t buf3[NTRUPLUS_SSBYTES + NTRUPLUS_N / 4];
	} u;

	int8_t fail;

	poly a;  /* [OPT8] Only 1 poly! */

	/* Step 1-2: c * f -> a (both from bytes) */
	poly_basemul_2tobytes(&a, ct, sk);

	/* Step 3: INVNTT + crepmod3 */
	poly_invntt_crepmod3(&a, &a);

	/* Step 4: compress m1 FIRST, then NTT in-place */
	poly_compress_ternary(u.m1_ternary, &a);
	poly_ntt_ternary(&a, &a);  /* [OPT8] in-place NTT! a = NTT(m1) */

	/* Steps 5-6 fused: (c - NTT(m1)) * hinv, overwriting a block by block */
	poly_sub_basemul_tobytes(&a, ct, &a, sk + NTRUPLUS_POLYBYTES);

	/* [OPT9] a = r2. Compute hash_g(r2) ONCE; reuse for sotp_inv AND verify
	 * (was computed twice into r2_hash and u.hashbuf — identical result). */
	poly_tobytes_hash_g(r2_hash, &a);

	/* Step 7: decode msg using m1 from compressed bytes directly */
	fail = poly_sotp_inv_compressed(msg, u.m1_ternary, r2_hash);

	/* Step 8: reconstruct hash key */
	for (int i = 0; i < NTRUPLUS_SYMBYTES; i++)
		msg[i + NTRUPLUS_N / 8] = sk[i + 2 * NTRUPLUS_POLYBYTES];

	/* Step 9: re-derive r1 */
	for (int i = 0; i < NTRUPLUS_SSBYTES + NTRUPLUS_N / 4; i++)
		u.buf3[i] = 0;

	hash_h_kem(u.buf3, msg);

	poly_cbd1(&a, u.buf3 + NTRUPLUS_SSBYTES);
	poly_ntt_ternary(&a, &a);

	/* [OPT8] Hash-based verify: compare hash_g(r1) with hash_g(r2) */
	{
		uint8_t r1_hash[NTRUPLUS_N / 4];
		poly_tobytes_hash_g(r1_hash, &a);
		fail |= verify(r1_hash, r2_hash, NTRUPLUS_N / 4);
	}

	for(int i = 0; i < NTRUPLUS_SSBYTES; i++)
		ss[i] = u.buf3[i] & ~(-fail);

	return fail;
}
