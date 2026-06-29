// SPDX-License-Identifier: MIT

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "symmetric.h"
#include "fips202.h"  /* [OPT3] for shake256_inc_* in poly_tobytes_hash_g */

/*************************************************
* Name:        load16_littleendian
**************************************************/
static uint16_t load16_littleendian(const uint8_t x[2])
{
	uint16_t r;
	r  = (uint32_t)x[0];
	r |= (uint32_t)x[1] << 8;;
	return r;
}

/*************************************************
* Name:        crepmod3
**************************************************/
static int16_t crepmod3(int16_t a)
{
	a += (a >> 15) & NTRUPLUS_Q;
	a -= (NTRUPLUS_Q+1)/2;
	a += (a >> 15) & NTRUPLUS_Q;
	a -= (NTRUPLUS_Q-1)/2;

	a  = (a >> 8) + (a & 255);
	a  = (a >> 4) + (a & 15);
	a  = (a >> 2) + (a & 3);
	a  = (a >> 2) + (a & 3);
	a -= 3;
	a += ((a + 1) >> 15) & 3;
	return a;
}

/*************************************************
* Name:        poly_tobytes
**************************************************/
void poly_tobytes(uint8_t r[NTRUPLUS_POLYBYTES], const poly *a)
{
	int16_t t[4];

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 18; j++)
		{
			t[0]  = a->coeffs[64*j + i];
			t[0] += (t[0] >> 15) & NTRUPLUS_Q;
			t[1]  = a->coeffs[64*j + i + 16];
			t[1] += (t[1] >> 15) & NTRUPLUS_Q;
			t[2]  = a->coeffs[64*j + i + 32];
			t[2] += (t[2] >> 15) & NTRUPLUS_Q;
			t[3]  = a->coeffs[64*j + i + 48];
			t[3] += (t[3] >> 15) & NTRUPLUS_Q;

			r[96*j + 2*i +  0] = t[0];
			r[96*j + 2*i +  1] = (t[0] >> 8) | (t[1] << 4);
			r[96*j + 2*i + 32] = (t[1] >> 4);
			r[96*j + 2*i + 33] = t[2];
			r[96*j + 2*i + 64] = (t[2] >> 8) | (t[3] << 4);
			r[96*j + 2*i + 65] = (t[3] >> 4);
		}
	}
}

/*************************************************
* Name:        poly_frombytes
**************************************************/
void poly_frombytes(poly *r, const uint8_t a[NTRUPLUS_POLYBYTES])
{
	unsigned char t[6];

	for(int i = 0; i < 16; i++)
	{
		for(int j = 0; j < 18; j++)
		{
			t[0] = a[96*j + 2*i];
			t[1] = a[96*j + 2*i + 1];
			t[2] = a[96*j + 2*i + 32];
			t[3] = a[96*j + 2*i + 33];
			t[4] = a[96*j + 2*i + 64];
			t[5] = a[96*j + 2*i + 65];

			r->coeffs[64*j + i +  0] = t[0]      | ((int16_t)t[1] & 0xf) << 8;
			r->coeffs[64*j + i + 16] = t[1] >> 4 | ((int16_t)t[2]      ) << 4;
			r->coeffs[64*j + i + 32] = t[3]      | ((int16_t)t[4] & 0xf) << 8;
			r->coeffs[64*j + i + 48] = t[4] >> 4 | ((int16_t)t[5]      ) << 4;
		}
	}
}

/*************************************************
* Name:        poly_cbd1
**************************************************/
void poly_cbd1(poly *r, const unsigned char buf[NTRUPLUS_N/4])
{
	uint16_t t1, t2;

	for(int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			t1 = load16_littleendian(buf + 32*i + 2*j);
			t2 = load16_littleendian(buf + 32*i + 2*j + 144);

			for(int k = 0; k < 16; k++)
			{
				r->coeffs[256*i + 16*k + j] = (t1 & 0x1) - (t2 & 0x1);

				t1 >>= 1;
				t2 >>= 1;
			}
		}
	}

	for (int j = 0; j < 8; j++)
	{
		t1 = load16_littleendian(buf + 128 + 2*j);
		t2 = load16_littleendian(buf + 128 + 2*j + 144);

		for(int k = 0; k < 16; k++)
		{
			r->coeffs[1024 + 8*k + j] = (t1 & 0x1) - (t2 & 0x1);

			t1 >>= 1;
			t2 >>= 1;
		}
	}
}

/*************************************************
* [OPT9] poly_cbd1_from_seed
*
* Streaming CBD1: generates poly directly from 32-byte seed.
* Squeezes SHAKE256 in two 144-byte halves.
**************************************************/
void poly_cbd1_from_seed(poly *r, const uint8_t seed[32])
{
	shake256incctx ctx;
	uint8_t half[144];
	uint16_t t;

	shake256_inc_init(&ctx);
	shake256_inc_absorb(&ctx, seed, 32);
	shake256_inc_finalize(&ctx);

	/* Phase 1: first 144 bytes -> positive */
	shake256_inc_squeeze(half, 144, &ctx);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 16; j++) {
			t = load16_littleendian(half + 32*i + 2*j);
			for (int k = 0; k < 16; k++) {
				r->coeffs[256*i + 16*k + j] = (t & 0x1);
				t >>= 1;
			}
		}
	for (int j = 0; j < 8; j++) {
		t = load16_littleendian(half + 128 + 2*j);
		for (int k = 0; k < 16; k++) {
			r->coeffs[1024 + 8*k + j] = (t & 0x1);
			t >>= 1;
		}
	}

	/* Phase 2: next 144 bytes -> subtract negative */
	shake256_inc_squeeze(half, 144, &ctx);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 16; j++) {
			t = load16_littleendian(half + 32*i + 2*j);
			for (int k = 0; k < 16; k++) {
				r->coeffs[256*i + 16*k + j] -= (t & 0x1);
				t >>= 1;
			}
		}
	for (int j = 0; j < 8; j++) {
		t = load16_littleendian(half + 128 + 2*j);
		for (int k = 0; k < 16; k++) {
			r->coeffs[1024 + 8*k + j] -= (t & 0x1);
			t >>= 1;
		}
	}

	shake256_inc_ctx_release(&ctx);
}

/*************************************************
* Name:        poly_sotp
**************************************************/
void poly_sotp(poly *r, const uint8_t *msg, const uint8_t *buf)
{
    uint8_t tmp[NTRUPLUS_N/4];

    for(int i = 0; i < NTRUPLUS_N/8; i++)
    {
         tmp[i] = buf[i]^msg[i];
    }

    for(int i = NTRUPLUS_N/8; i < NTRUPLUS_N/4; i++)
    {
         tmp[i] = buf[i];
    }

	poly_cbd1(r, tmp);
}

/*************************************************
* Name:        poly_sotp_inv
**************************************************/
int poly_sotp_inv(unsigned char *msg, const poly *a, const unsigned char *buf)
{
	uint16_t t1, t2, t3, t4;
	uint32_t r = 0;

	for(int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			t1 = load16_littleendian(buf + 32*i + 2*j);
			t2 = load16_littleendian(buf + 32*i + 2*j + 144);
			t3 = 0;

			for(int k = 0; k < 16; k++)
			{
				t4 = t2 & 0x1;
				t4 += a->coeffs[256*i + 16*k + j];
				r |= t4;
				t4 = (t4^t1) & 0x1;
				t3 ^= t4 << k;

				t1 >>= 1;
				t2 >>= 1;
			}

			msg[32*i + 2*j   ] = t3;
			msg[32*i + 2*j + 1] = t3 >> 8;
		}
	}

	for (int j = 0; j < 8; j++)
	{
		t1 = load16_littleendian(buf + 128 + 2*j);
		t2 = load16_littleendian(buf + 128 + 2*j + 144);
		t3 = 0;

		for(int k = 0; k < 16; k++)
		{
			t4 = t2 & 0x1;
			t4 += a->coeffs[1024 + 8*k + j];
			r |= t4;
			t4 = (t4^t1) & 0x1;
			t3 ^= t4 << k;

			t1 >>= 1;
			t2 >>= 1;
		}

		msg[128 + 2*j   ] = t3;
		msg[128 + 2*j + 1] = t3 >> 8;
	}

	r = r >> 1;
	r = (-(uint64_t)r) >> 63;

	return r;
}

/*************************************************
* Name:        poly_ntt
**************************************************/
void poly_ntt(poly *r, const poly *a)
{
	ntt(r->coeffs, a->coeffs);
}

/* [OPT-NTT] NTT for ternary input: skips Barrett reduction */
void poly_ntt_ternary(poly *r, const poly *a)
{
	ntt(r->coeffs, a->coeffs);          /* [PORT] ternary->full ntt */
}

/*************************************************
* Name:        poly_invntt
**************************************************/
void poly_invntt(poly *r, const poly *a)
{
	invntt(r->coeffs, a->coeffs);
}

/*************************************************
* Name:        poly_baseinv
**************************************************/
int poly_baseinv(poly *r, const poly *a)
{
	int result = 0;

	for(int i = 0; i < NTRUPLUS_N/8; ++i)
	{
		result = baseinv(r->coeffs + 8*i, a->coeffs + 8*i, zetas[144 + i]);
		if(result) return 1;
		result = baseinv(r->coeffs + 8*i + 4, a->coeffs + 8*i + 4, -zetas[144 + i]);
		if(result) return 1;
	 }

	return 0;
}

/*************************************************
* Name:        poly_basemul
**************************************************/
void poly_basemul(poly *r, const poly *a, const poly *b)
{
	for(int i = 0; i < NTRUPLUS_N/8; ++i)
	{
		basemul(r->coeffs + 8*i, a->coeffs + 8*i, b->coeffs + 8*i, zetas[144 + i]);
		basemul(r->coeffs + 8*i + 4, a->coeffs + 8*i + 4, b->coeffs + 8*i + 4, -zetas[144 + i]);
	}
}

/*************************************************
* Name:        poly_basemul_add
**************************************************/
void poly_basemul_add(poly *r, const poly *a, const poly *b, const poly *c)
{
	for(int i = 0; i < NTRUPLUS_N/8; ++i)
	{
		basemul_add(r->coeffs + 8*i, a->coeffs + 8*i, b->coeffs + 8*i, c->coeffs + 8*i, zetas[144 + i]);
		basemul_add(r->coeffs + 8*i + 4, a->coeffs + 8*i + 4, b->coeffs + 8*i + 4, c->coeffs + 8*i + 4, -zetas[144 + i]);
	}
}

/*************************************************
* Name:        poly_sub
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
	for(int i = 0; i < NTRUPLUS_N; ++i)
		r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

/*************************************************
* Name:        poly_triple
**************************************************/
void poly_triple(poly *r, const poly *a)
{
	for(int i = 0; i < NTRUPLUS_N; ++i)
		r->coeffs[i] = 3*a->coeffs[i];
}

/*************************************************
* Name:        poly_crepmod3
**************************************************/
void poly_crepmod3(poly *r, const poly *a)
{
	for(int i = 0; i < NTRUPLUS_N; i++)
		r->coeffs[i] = crepmod3(a->coeffs[i]);
}

/* [OPT6] wide_crepmod3: centered mod 3 for arbitrary int16_t.
 * Works for any Q where Q != 0 (mod 3), without Q-dependent normalization.
 * Uses uint16_t reinterpretation (65535 = 0 mod 3) + base-256/16/4 digit sum. */
static int16_t wide_crepmod3(int16_t a)
{
	uint16_t u = (uint16_t)(a - ((a >> 15) & 1));

	uint16_t s = (u >> 8) + (u & 0xFF);
	s = (s >> 4) + (s & 0xF);
	s = (s >> 2) + (s & 0x3);
	s = (s >> 2) + (s & 0x3);
	s = (s >> 2) + (s & 0x3);

	int16_t r2 = (int16_t)s - 3;
	r2 += ((r2 + 1) >> 15) & 3;
	return r2;
}

/* [OPT6] Fused invntt_lazy + crepmod3. */
void poly_invntt_crepmod3(poly *r, const poly *a)
{
	invntt(r->coeffs, a->coeffs);        /* [PORT] lazy->full invntt */

	for(int i = 0; i < NTRUPLUS_N; i++)
		r->coeffs[i] = wide_crepmod3(r->coeffs[i]);
}

/*************************************************
* [OPT7] poly_compress_ternary
*
* Compress ternary polynomial ({-1,0,1} coefficients) to 2 bits/coeff.
* Encoding: -1 -> 0b11, 0 -> 0b00, 1 -> 0b01
* Output: N/4 = 288 bytes (vs 2304 bytes for int16_t[N])
**************************************************/
void poly_compress_ternary(uint8_t r[NTRUPLUS_N/4], const poly *a)
{
	for (int i = 0; i < NTRUPLUS_N / 4; i++)
	{
		uint8_t byte = 0;
		for (int j = 0; j < 4; j++)
		{
			byte |= ((uint8_t)(a->coeffs[4*i + j] & 0x3)) << (2*j);
		}
		r[i] = byte;
	}
}

/*************************************************
* [OPT7] poly_decompress_ternary
*
* Decompress ternary polynomial from 2-bit packed form.
* Decoding: 0b11(3) -> -1, 0b00(0) -> 0, 0b01(1) -> 1
**************************************************/
void poly_decompress_ternary(poly *r, const uint8_t a[NTRUPLUS_N/4])
{
	for (int i = 0; i < NTRUPLUS_N / 4; i++)
	{
		uint8_t byte = a[i];
		for (int j = 0; j < 4; j++)
		{
			uint8_t bits = (byte >> (2*j)) & 0x3;
			int16_t c = (int16_t)bits;
			if (c == 3) c = -1;
			r->coeffs[4*i + j] = c;
		}
	}
}

/*************************************************
* [OPT8] poly_pack_seq
*
* Sequential 12-bit packing: 2 coefficients per 3 bytes.
* Output size: N * 12 / 8 = 1728 bytes = POLYBYTES.
**************************************************/
void poly_pack_seq(uint8_t r[NTRUPLUS_POLYBYTES], const poly *a)
{
	for (int i = 0; i < NTRUPLUS_N / 2; i++)
	{
		int16_t t0 = a->coeffs[2*i];
		int16_t t1 = a->coeffs[2*i + 1];
		t0 += (t0 >> 15) & NTRUPLUS_Q;
		t1 += (t1 >> 15) & NTRUPLUS_Q;
		r[3*i + 0] = (uint8_t)t0;
		r[3*i + 1] = (uint8_t)((t0 >> 8) | (t1 << 4));
		r[3*i + 2] = (uint8_t)(t1 >> 4);
	}
}

/*************************************************
* [OPT8] Unpack 8 sequential coefficients from packed form.
* 8 coefficients * 12 bits = 96 bits = 12 bytes per chunk.
**************************************************/
static void unpack_8_coeffs(int16_t r[8], const uint8_t *p)
{
	for (int i = 0; i < 4; i++)
	{
		r[2*i]     =  (int16_t)p[3*i]          | (((int16_t)p[3*i+1] & 0xF) << 8);
		r[2*i + 1] = ((int16_t)p[3*i+1] >> 4)  | ( (int16_t)p[3*i+2]        << 4);
	}
}

/*************************************************
* [OPT8] poly_basemul_add_packed
*
* basemul_add where the second operand (b) is in sequential
* 12-bit packed form. Unpacks 8 coefficients at a time.
* r = a * b_packed + c  (in NTT domain, component-wise)
**************************************************/
void poly_basemul_add_packed(poly *r, const poly *a, const uint8_t *b_packed, const poly *c)
{
	int16_t b_chunk[8];

	for (int i = 0; i < NTRUPLUS_N / 8; ++i)
	{
		unpack_8_coeffs(b_chunk, b_packed + 12*i);

		basemul_add(r->coeffs + 8*i,     a->coeffs + 8*i,
		            b_chunk,              c->coeffs + 8*i,      zetas[144 + i]);
		basemul_add(r->coeffs + 8*i + 4, a->coeffs + 8*i + 4,
		            b_chunk + 4,          c->coeffs + 8*i + 4, -zetas[144 + i]);
	}
}

/* [OPT5] Unpack 8 coefficients from tobytes interleaved format by basemul block.
 * N=1152: 18 groups of 64, no tail. stride=8 always within one sub-group. */
static void unpack_block_from_tobytes(int16_t r[8], const uint8_t *a, int block_idx)
{
	int start = 8 * block_idx;
	int j = start >> 6;
	int sub = (start >> 4) & 3;
	int i0 = start & 15;
	const uint8_t *p = a + 96*j + 2*i0;

	switch (sub) {
	case 0:
		for (int k = 0; k < 8; k++)
			r[k] =  p[2*k]       | (((int16_t)p[2*k+ 1] & 0xF) << 8);
		break;
	case 1:
		for (int k = 0; k < 8; k++)
			r[k] = (p[2*k+ 1]>>4)| ( (int16_t)p[2*k+32]        << 4);
		break;
	case 2:
		for (int k = 0; k < 8; k++)
			r[k] =  p[2*k+33]    | (((int16_t)p[2*k+64] & 0xF) << 8);
		break;
	case 3:
		for (int k = 0; k < 8; k++)
			r[k] = (p[2*k+64]>>4)| ( (int16_t)p[2*k+65]        << 4);
		break;
	}
}

/* [OPT7] poly_basemul_packed_tobytes: one operand from tobytes format.
 * Supports r == a aliasing (basemul is aliasing-safe via OPT1). */
void poly_basemul_packed_tobytes(poly *r, const poly *a, const uint8_t *b_tobytes)
{
	int16_t b_chunk[8];

	for (int i = 0; i < NTRUPLUS_N / 8; ++i)
	{
		unpack_block_from_tobytes(b_chunk, b_tobytes, i);

		basemul(r->coeffs + 8*i,     a->coeffs + 8*i,
		        b_chunk,              zetas[144 + i]);
		basemul(r->coeffs + 8*i + 4, a->coeffs + 8*i + 4,
		        b_chunk + 4,         -zetas[144 + i]);
	}
}

/* [OPT5] poly_basemul_add_2packed: both operands from packed bytes.
 * Eliminates poly b in encaps, saving N*2 bytes of stack. */
void poly_basemul_add_2packed(poly *r, const uint8_t *a_tobytes,
                              const uint8_t *b_packseq, const poly *c)
{
	int16_t a_chunk[8], b_chunk[8];

	for (int i = 0; i < NTRUPLUS_N / 8; ++i)
	{
		unpack_block_from_tobytes(a_chunk, a_tobytes, i);
		unpack_8_coeffs(b_chunk, b_packseq + 12*i);

		basemul_add(r->coeffs + 8*i,     a_chunk,
		            b_chunk,              c->coeffs + 8*i,      zetas[144 + i]);
		basemul_add(r->coeffs + 8*i + 4, a_chunk + 4,
		            b_chunk + 4,          c->coeffs + 8*i + 4, -zetas[144 + i]);
	}
}

/*************************************************
* [OPT3] poly_tobytes_hash_g
*
* Fused serialization + hash_g: serializes polynomial
* block-by-block into SHAKE256 incremental absorb.
* Uses j<18 for 1152 (N=1152, 1152/64=18 blocks of 96 bytes).
**************************************************/
void poly_tobytes_hash_g(uint8_t *hash_out, const poly *a)
{
	shake256incctx ctx;
	uint8_t domain = 0x1;
	uint8_t chunk[96];

	shake256_inc_init(&ctx);
	shake256_inc_absorb(&ctx, &domain, 1);

	for (int j = 0; j < 18; j++)
	{
		int16_t t[4];
		for (int i = 0; i < 16; i++)
		{
			t[0]  = a->coeffs[64*j + i];
			t[0] += (t[0] >> 15) & NTRUPLUS_Q;
			t[1]  = a->coeffs[64*j + i + 16];
			t[1] += (t[1] >> 15) & NTRUPLUS_Q;
			t[2]  = a->coeffs[64*j + i + 32];
			t[2] += (t[2] >> 15) & NTRUPLUS_Q;
			t[3]  = a->coeffs[64*j + i + 48];
			t[3] += (t[3] >> 15) & NTRUPLUS_Q;

			chunk[2*i +  0] = t[0];
			chunk[2*i +  1] = (t[0] >> 8) | (t[1] << 4);
			chunk[2*i + 32] = (t[1] >> 4);
			chunk[2*i + 33] = t[2];
			chunk[2*i + 64] = (t[2] >> 8) | (t[3] << 4);
			chunk[2*i + 65] = (t[3] >> 4);
		}
		shake256_inc_absorb(&ctx, chunk, 96);
	}

	shake256_inc_finalize(&ctx);
	shake256_inc_squeeze(hash_out, NTRUPLUS_N / 4, &ctx);
	shake256_inc_ctx_release(&ctx);
}

/*************************************************
* [OPT3] poly_verify
*
* Constant-time comparison of two polynomials.
**************************************************/
int poly_verify(const poly *a, const poly *b)
{
	uint16_t r = 0;

	for (int i = 0; i < NTRUPLUS_N; i++)
	{
		int16_t ca = a->coeffs[i];
		int16_t cb = b->coeffs[i];
		ca += (ca >> 15) & NTRUPLUS_Q;
		cb += (cb >> 15) & NTRUPLUS_Q;
		r |= (uint16_t)((uint16_t)ca ^ (uint16_t)cb);
	}

	return (-(uint64_t)r) >> 63;
}

/* [OPT8] poly_basemul_2tobytes */
void poly_basemul_2tobytes(poly *r, const uint8_t *a_tobytes, const uint8_t *b_tobytes)
{
	int16_t a_chunk[8], b_chunk[8];
	for (int i = 0; i < NTRUPLUS_N / 8; ++i)
	{
		unpack_block_from_tobytes(a_chunk, a_tobytes, i);
		unpack_block_from_tobytes(b_chunk, b_tobytes, i);
		basemul(r->coeffs + 8*i,     a_chunk,     b_chunk,      zetas[144 + i]);
		basemul(r->coeffs + 8*i + 4, a_chunk + 4, b_chunk + 4, -zetas[144 + i]);
	}
}

/* [OPT8] poly_sub_basemul_tobytes */
void poly_sub_basemul_tobytes(poly *r, const uint8_t *ct_tobytes,
                               const poly *a, const uint8_t *hinv_tobytes)
{
	int16_t ct_chunk[8], hinv_chunk[8], diff[8];
	for (int i = 0; i < NTRUPLUS_N / 8; ++i)
	{
		unpack_block_from_tobytes(ct_chunk, ct_tobytes, i);
		unpack_block_from_tobytes(hinv_chunk, hinv_tobytes, i);
		for (int k = 0; k < 8; k++)
			diff[k] = ct_chunk[k] - a->coeffs[8*i + k];
		basemul(r->coeffs + 8*i,     diff,     hinv_chunk,      zetas[144 + i]);
		basemul(r->coeffs + 8*i + 4, diff + 4, hinv_chunk + 4, -zetas[144 + i]);
	}
}

/* [OPT8] poly_sotp_inv_compressed */
static int16_t decompress_one_ternary(const uint8_t *m1t, int idx)
{
	uint8_t bits = (m1t[idx >> 2] >> (2 * (idx & 3))) & 0x3;
	return (bits == 3) ? -1 : (int16_t)bits;
}

int poly_sotp_inv_compressed(unsigned char *msg, const uint8_t *m1t, const unsigned char *buf)
{
	uint16_t t1, t2, t3, t4;
	uint32_t r = 0;

	for(int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			t1 = load16_littleendian(buf + 32*i + 2*j);
			t2 = load16_littleendian(buf + 32*i + 2*j + 144);
			t3 = 0;

			for(int k = 0; k < 16; k++)
			{
				t4 = t2 & 0x1;
				t4 += decompress_one_ternary(m1t, 256*i + 16*k + j);
				r |= t4;
				t4 = (t4^t1) & 0x1;
				t3 ^= t4 << k;
				t1 >>= 1;
				t2 >>= 1;
			}
			msg[32*i + 2*j   ] = t3;
			msg[32*i + 2*j + 1] = t3 >> 8;
		}
	}

	for (int j = 0; j < 8; j++)
	{
		t1 = load16_littleendian(buf + 128 + 2*j);
		t2 = load16_littleendian(buf + 128 + 2*j + 144);
		t3 = 0;

		for(int k = 0; k < 16; k++)
		{
			t4 = t2 & 0x1;
			t4 += decompress_one_ternary(m1t, 1024 + 8*k + j);
			r |= t4;
			t4 = (t4^t1) & 0x1;
			t3 ^= t4 << k;
			t1 >>= 1;
			t2 >>= 1;
		}
		msg[128 + 2*j   ] = t3;
		msg[128 + 2*j + 1] = t3 >> 8;
	}

	r = r >> 1;
	r = (-(uint64_t)r) >> 63;
	return r;
}
