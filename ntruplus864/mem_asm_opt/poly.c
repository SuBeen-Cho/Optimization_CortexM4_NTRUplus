// SPDX-License-Identifier: MIT

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "symmetric.h"
#include "fips202.h"

static uint16_t load16_littleendian(const uint8_t x[2])
{
	uint16_t r;
	r  = (uint32_t)x[0];
	r |= (uint32_t)x[1] << 8;
	return r;
}

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

void poly_tobytes(uint8_t r[NTRUPLUS_POLYBYTES], const poly *a)
{
	int16_t t[4];

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 13; j++)
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

	for (int i = 0; i < 8; i++)
	{
		t[0]  = a->coeffs[832 + i];
		t[0] += (t[0] >> 15) & NTRUPLUS_Q;
		t[1]  = a->coeffs[832 + i + 8];
		t[1] += (t[1] >> 15) & NTRUPLUS_Q;
		t[2]  = a->coeffs[832 + i + 16];
		t[2] += (t[2] >> 15) & NTRUPLUS_Q;
		t[3]  = a->coeffs[832 + i + 24];
		t[3] += (t[3] >> 15) & NTRUPLUS_Q;

		r[1248 + 2*i +  0] = t[0];
		r[1248 + 2*i +  1] = (t[0] >> 8) | (t[1] << 4);
		r[1248 + 2*i + 16] = (t[1] >> 4);
		r[1248 + 2*i + 17] = t[2];
		r[1248 + 2*i + 32] = (t[2] >> 8) | (t[3] << 4);
		r[1248 + 2*i + 33] = (t[3] >> 4);
	}
}

void poly_frombytes(poly *r, const uint8_t a[NTRUPLUS_POLYBYTES])
{
	unsigned char t[6];

	for(int i = 0; i < 16; i++)
	{
		for(int j = 0; j < 13; j++)
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

	for(int i = 0; i < 8; i++)
	{
		t[0] = a[1248 + 2*i];
		t[1] = a[1248 + 2*i + 1];
		t[2] = a[1248 + 2*i + 16];
		t[3] = a[1248 + 2*i + 17];
		t[4] = a[1248 + 2*i + 32];
		t[5] = a[1248 + 2*i + 33];

		r->coeffs[832 + i +  0] = t[0]      | ((int16_t)t[1] & 0xf) << 8;
		r->coeffs[832 + i +  8] = t[1] >> 4 | ((int16_t)t[2]      ) << 4;
		r->coeffs[832 + i + 16] = t[3]      | ((int16_t)t[4] & 0xf) << 8;
		r->coeffs[832 + i + 24] = t[4] >> 4 | ((int16_t)t[5]      ) << 4;
	}
}

void poly_cbd1(poly *r, const unsigned char buf[NTRUPLUS_N/4])
{
	uint16_t t1, t2;

	for(int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			t1 = load16_littleendian(buf + 32*i + 2*j);
			t2 = load16_littleendian(buf + 32*i + 2*j + 108);

			for(int k = 0; k < 16; k++)
			{
				r->coeffs[256*i + 16*k + j] = (t1 & 0x1) - (t2 & 0x1);

				t1 >>= 1;
				t2 >>= 1;
			}
		}
	}

	for (int j = 0; j < 4; j++)
	{
		t1 = load16_littleendian(buf + 96 + 2*j);
		t2 = load16_littleendian(buf + 96 + 2*j + 108);

		for(int k = 0; k < 16; k++)
		{
			r->coeffs[768 + 4*k + j] = (t1 & 0x1) - (t2 & 0x1);

			t1 >>= 1;
			t2 >>= 1;
		}
	}

	for (int j = 0; j < 2; j++)
	{
		t1 = load16_littleendian(buf + 104 + 2*j);
		t2 = load16_littleendian(buf + 104 + 2*j + 108);

		for(int k = 0; k < 16; k++)
		{
			r->coeffs[832 + 2*k + j] = (t1 & 0x1) - (t2 & 0x1);

			t1 >>= 1;
			t2 >>= 1;
		}
	}
}

/*
 * [OPT9] poly_cbd1_from_seed
 * Streaming CBD1: generates poly directly from 32-byte seed.
 * N=864, N/8=108: squeezes 108 bytes per half.
 */
void poly_cbd1_from_seed(poly *r, const uint8_t seed[32])
{
	shake256incctx ctx;
	uint8_t half[108];  /* N/8 = 108 */
	uint16_t t;

	shake256_inc_init(&ctx);
	shake256_inc_absorb(&ctx, seed, 32);
	shake256_inc_finalize(&ctx);

	/* Phase 1: first 108 bytes -> positive bits */
	shake256_inc_squeeze(half, 108, &ctx);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 16; j++) {
			t = load16_littleendian(half + 32*i + 2*j);
			for (int k = 0; k < 16; k++) {
				r->coeffs[256*i + 16*k + j] = (t & 0x1);
				t >>= 1;
			}
		}
	for (int j = 0; j < 4; j++) {
		t = load16_littleendian(half + 96 + 2*j);
		for (int k = 0; k < 16; k++) {
			r->coeffs[768 + 4*k + j] = (t & 0x1);
			t >>= 1;
		}
	}
	for (int j = 0; j < 2; j++) {
		t = load16_littleendian(half + 104 + 2*j);
		for (int k = 0; k < 16; k++) {
			r->coeffs[832 + 2*k + j] = (t & 0x1);
			t >>= 1;
		}
	}

	/* Phase 2: next 108 bytes -> subtract negative bits */
	shake256_inc_squeeze(half, 108, &ctx);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 16; j++) {
			t = load16_littleendian(half + 32*i + 2*j);
			for (int k = 0; k < 16; k++) {
				r->coeffs[256*i + 16*k + j] -= (t & 0x1);
				t >>= 1;
			}
		}
	for (int j = 0; j < 4; j++) {
		t = load16_littleendian(half + 96 + 2*j);
		for (int k = 0; k < 16; k++) {
			r->coeffs[768 + 4*k + j] -= (t & 0x1);
			t >>= 1;
		}
	}
	for (int j = 0; j < 2; j++) {
		t = load16_littleendian(half + 104 + 2*j);
		for (int k = 0; k < 16; k++) {
			r->coeffs[832 + 2*k + j] -= (t & 0x1);
			t >>= 1;
		}
	}

	shake256_inc_ctx_release(&ctx);
}

void poly_sotp(poly *r, const uint8_t *msg, const uint8_t *buf)
{
	uint8_t tmp[NTRUPLUS_N/4];

	for(int i = 0; i < NTRUPLUS_N/8; i++)
		tmp[i] = buf[i]^msg[i];

	for(int i = NTRUPLUS_N/8; i < NTRUPLUS_N/4; i++)
		tmp[i] = buf[i];

	poly_cbd1(r, tmp);
}

int poly_sotp_inv(unsigned char *msg, const poly *a, const unsigned char *buf)
{
	uint32_t t1, t2, t3, t4;
	uint32_t r = 0;

	for(int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			t1 = load16_littleendian(buf + 32*i + 2*j);
			t2 = load16_littleendian(buf + 32*i + 2*j + 108);
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

			msg[32*i + 2*j    ] = t3;
			msg[32*i + 2*j + 1] = t3 >> 8;
		}
	}

	for (int j = 0; j < 4; j++)
	{
		t1 = load16_littleendian(buf + 96 + 2*j);
		t2 = load16_littleendian(buf + 96 + 2*j + 108);
		t3 = 0;

		for(int k = 0; k < 16; k++)
		{
			t4 = t2 & 0x1;
			t4 += a->coeffs[768 + 4*k + j];
			r |= t4;
			t4 = (t4^t1) & 0x1;
			t3 ^= t4 << k;

			t1 >>= 1;
			t2 >>= 1;
		}

		msg[96 + 2*j    ] = t3;
		msg[96 + 2*j + 1] = t3 >> 8;
	}

	for (int j = 0; j < 2; j++)
	{
		t1 = load16_littleendian(buf + 104 + 2*j);
		t2 = load16_littleendian(buf + 104 + 2*j + 108);
		t3 = 0;

		for(int k = 0; k < 16; k++)
		{
			t4 = t2 & 0x1;
			t4 += a->coeffs[832 + 2*k + j];
			r |= t4;
			t4 = (t4^t1) & 0x1;
			t3 ^= t4 << k;

			t1 >>= 1;
			t2 >>= 1;
		}

		msg[104 + 2*j    ] = t3;
		msg[104 + 2*j + 1] = t3 >> 8;
	}

	r = r >> 1;
	r = (-(uint64_t)r) >> 63;

	return r;
}

void poly_ntt(poly *r, const poly *a)
{
	ntt(r->coeffs, a->coeffs);
}

/* [OPT-NTT] NTT for ternary input: skips Barrett reduction */
void poly_ntt_ternary(poly *r, const poly *a)
{
	ntt(r->coeffs, a->coeffs);          /* [PORT] ternary->full ntt */
}

void poly_invntt(poly *r, const poly *a)
{
	invntt(r->coeffs, a->coeffs);
}

int poly_baseinv(poly *r, const poly *a)
{
	int result = 0;

	for(int i = 0; i < NTRUPLUS_N/6; ++i)
	{
		result = baseinv(r->coeffs+6*i, a->coeffs+6*i, zetas[144+i]);
		if(result) return 1;
		result = baseinv(r->coeffs+6*i+3, a->coeffs+6*i+3, -zetas[144+i]);
		if(result) return 1;
	}

	return 0;
}

void poly_basemul(poly *r, const poly *a, const poly *b)
{
	for(int i = 0; i < NTRUPLUS_N/6; i++)
	{
		basemul(r->coeffs+6*i,   a->coeffs+6*i,   b->coeffs+6*i,    zetas[144+i]);
		basemul(r->coeffs+6*i+3, a->coeffs+6*i+3, b->coeffs+6*i+3, -zetas[144+i]);
	}
}

void poly_basemul_add(poly *r, const poly *a, const poly *b, const poly *c)
{
	for(int i = 0; i < NTRUPLUS_N/6; i++)
	{
		basemul_add(r->coeffs+6*i,   a->coeffs+6*i,   b->coeffs+6*i,   c->coeffs+6*i,    zetas[144+i]);
		basemul_add(r->coeffs+6*i+3, a->coeffs+6*i+3, b->coeffs+6*i+3, c->coeffs+6*i+3, -zetas[144+i]);
	}
}

void poly_sub(poly *r, const poly *a, const poly *b)
{
	for(int i = 0; i < NTRUPLUS_N; ++i)
		r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

void poly_triple(poly *r, const poly *a)
{
	for(int i = 0; i < NTRUPLUS_N; ++i)
		r->coeffs[i] = 3*a->coeffs[i];
}

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

/*
 * [OPT7] poly_compress_ternary
 * N/4 = 216 bytes for 864
 */
void poly_compress_ternary(uint8_t r[NTRUPLUS_N/4], const poly *a)
{
	for (int i = 0; i < NTRUPLUS_N / 4; i++)
	{
		uint8_t byte = 0;
		for (int j = 0; j < 4; j++)
			byte |= ((uint8_t)(a->coeffs[4*i + j] & 0x3)) << (2*j);
		r[i] = byte;
	}
}

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

/*
 * [OPT8] poly_pack_seq
 * Sequential 12-bit packing: 2 coefficients per 3 bytes.
 * N=864: POLYBYTES=1296 = 864*12/8
 */
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

/*
 * [OPT8] Unpack 6 sequential coefficients (9 bytes).
 * 6 coefficients * 12 bits = 72 bits = 9 bytes
 */
static void unpack_6_coeffs(int16_t r[6], const uint8_t *p)
{
	for (int i = 0; i < 3; i++)
	{
		r[2*i]     =  (int16_t)p[3*i]          | (((int16_t)p[3*i+1] & 0xF) << 8);
		r[2*i + 1] = ((int16_t)p[3*i+1] >> 4)  | ( (int16_t)p[3*i+2]        << 4);
	}
}

/*
 * [OPT8] poly_basemul_add_packed
 * For ntruplus864: 6-stride, N/6=144 blocks.
 * Each block: unpack 6 coefficients (9 bytes), call basemul_add on 3+3.
 */
void poly_basemul_add_packed(poly *r, const poly *a, const uint8_t *b_packed, const poly *c)
{
	int16_t b_chunk[6];

	for (int i = 0; i < NTRUPLUS_N / 6; ++i)
	{
		unpack_6_coeffs(b_chunk, b_packed + 9*i);

		basemul_add(r->coeffs + 6*i,   a->coeffs + 6*i,
		            b_chunk,            c->coeffs + 6*i,    zetas[144 + i]);
		basemul_add(r->coeffs + 6*i+3, a->coeffs + 6*i+3,
		            b_chunk + 3,        c->coeffs + 6*i+3, -zetas[144 + i]);
	}
}

/* [OPT5] Unpack one coefficient from tobytes interleaved format.
 * N=864: 13 groups of 64 (main) + 32 tail coefficients. */
static int16_t unpack_one_tobytes(const uint8_t *a, int idx)
{
	if (idx < 832)
	{
		int j   = idx >> 6;          /* idx / 64 */
		int sub = (idx >> 4) & 3;    /* (idx % 64) / 16 */
		int i   = idx & 15;         /* idx % 16 */
		const uint8_t *p = a + 96*j;

		switch (sub) {
		case 0: return  p[2*i]       | (((int16_t)p[2*i+ 1] & 0xF) << 8);
		case 1: return (p[2*i+ 1]>>4)| ( (int16_t)p[2*i+32]        << 4);
		case 2: return  p[2*i+33]    | (((int16_t)p[2*i+64] & 0xF) << 8);
		default:return (p[2*i+64]>>4)| ( (int16_t)p[2*i+65]        << 4);
		}
	}
	else
	{
		int tail = idx - 832;
		int sub  = tail >> 3;        /* tail / 8 */
		int i    = tail & 7;        /* tail % 8 */
		const uint8_t *p = a + 1248;

		switch (sub) {
		case 0: return  p[2*i]       | (((int16_t)p[2*i+ 1] & 0xF) << 8);
		case 1: return (p[2*i+ 1]>>4)| ( (int16_t)p[2*i+16]        << 4);
		case 2: return  p[2*i+17]    | (((int16_t)p[2*i+32] & 0xF) << 8);
		default:return (p[2*i+32]>>4)| ( (int16_t)p[2*i+33]        << 4);
		}
	}
}

/* [OPT7] poly_basemul_packed_tobytes: one operand from tobytes format.
 * For N=864: per-coefficient unpack (stride=6, blocks don't align with sub-groups).
 * Supports r == a aliasing (basemul is aliasing-safe via OPT1). */
void poly_basemul_packed_tobytes(poly *r, const poly *a, const uint8_t *b_tobytes)
{
	int16_t b_chunk[6];

	for (int i = 0; i < NTRUPLUS_N / 6; ++i)
	{
		for (int k = 0; k < 6; k++)
			b_chunk[k] = unpack_one_tobytes(b_tobytes, 6*i + k);

		basemul(r->coeffs + 6*i,     a->coeffs + 6*i,
		        b_chunk,              zetas[144 + i]);
		basemul(r->coeffs + 6*i + 3, a->coeffs + 6*i + 3,
		        b_chunk + 3,         -zetas[144 + i]);
	}
}

/* [OPT5] poly_basemul_add_2packed: both operands from packed bytes.
 * For N=864: stride=6, 144 blocks. Per-coefficient unpack for tobytes
 * because 6 doesn't divide 16 (blocks can span sub-groups). */
void poly_basemul_add_2packed(poly *r, const uint8_t *a_tobytes,
                              const uint8_t *b_packseq, const poly *c)
{
	int16_t a_chunk[6], b_chunk[6];

	for (int i = 0; i < NTRUPLUS_N / 6; ++i)
	{
		for (int k = 0; k < 6; k++)
			a_chunk[k] = unpack_one_tobytes(a_tobytes, 6*i + k);

		unpack_6_coeffs(b_chunk, b_packseq + 9*i);

		basemul_add(r->coeffs + 6*i,   a_chunk,
		            b_chunk,            c->coeffs + 6*i,    zetas[144 + i]);
		basemul_add(r->coeffs + 6*i+3, a_chunk + 3,
		            b_chunk + 3,        c->coeffs + 6*i+3, -zetas[144 + i]);
	}
}

/*
 * [OPT3] poly_tobytes_hash_g for ntruplus864 (N=864)
 * poly_tobytes structure: j<13 main (832 coeffs) + tail 32 coeffs.
 */
void poly_tobytes_hash_g(uint8_t *hash_out, const poly *a)
{
	shake256incctx ctx;
	uint8_t domain = 0x1;
	uint8_t chunk[96];

	shake256_inc_init(&ctx);
	shake256_inc_absorb(&ctx, &domain, 1);

	/* Main part: 13 blocks of 64 coefficients */
	for (int j = 0; j < 13; j++)
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

	/* Tail: 32 coefficients (indices 832..863), packed into 48 bytes */
	{
		uint8_t tail[48];
		int16_t t[4];
		for (int i = 0; i < 8; i++)
		{
			t[0] = a->coeffs[832 + i];
			t[0] += (t[0] >> 15) & NTRUPLUS_Q;
			t[1] = a->coeffs[832 + i + 8];
			t[1] += (t[1] >> 15) & NTRUPLUS_Q;
			t[2] = a->coeffs[832 + i + 16];
			t[2] += (t[2] >> 15) & NTRUPLUS_Q;
			t[3] = a->coeffs[832 + i + 24];
			t[3] += (t[3] >> 15) & NTRUPLUS_Q;
			tail[2*i +  0] = t[0];
			tail[2*i +  1] = (t[0] >> 8) | (t[1] << 4);
			tail[2*i + 16] = (t[1] >> 4);
			tail[2*i + 17] = t[2];
			tail[2*i + 32] = (t[2] >> 8) | (t[3] << 4);
			tail[2*i + 33] = (t[3] >> 4);
		}
		shake256_inc_absorb(&ctx, tail, 48);
	}

	shake256_inc_finalize(&ctx);
	shake256_inc_squeeze(hash_out, NTRUPLUS_N / 4, &ctx);
	shake256_inc_ctx_release(&ctx);
}

/*
 * [OPT3] poly_verify
 * Returns 0 if equal, 1 if different.
 */
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

/* [OPT8] poly_basemul_2tobytes: both operands from tobytes. */
void poly_basemul_2tobytes(poly *r, const uint8_t *a_tobytes, const uint8_t *b_tobytes)
{
	int16_t a_chunk[6], b_chunk[6];

	for (int i = 0; i < NTRUPLUS_N / 6; ++i)
	{
		for (int k = 0; k < 6; k++) {
			a_chunk[k] = unpack_one_tobytes(a_tobytes, 6*i + k);
			b_chunk[k] = unpack_one_tobytes(b_tobytes, 6*i + k);
		}

		basemul(r->coeffs + 6*i,     a_chunk,     b_chunk,      zetas[144 + i]);
		basemul(r->coeffs + 6*i + 3, a_chunk + 3, b_chunk + 3, -zetas[144 + i]);
	}
}

/* [OPT8] poly_sub_basemul_tobytes: fused (ct - a) * hinv from bytes. */
void poly_sub_basemul_tobytes(poly *r, const uint8_t *ct_tobytes,
                               const poly *a, const uint8_t *hinv_tobytes)
{
	int16_t ct_chunk[6], hinv_chunk[6], diff[6];

	for (int i = 0; i < NTRUPLUS_N / 6; ++i)
	{
		for (int k = 0; k < 6; k++) {
			ct_chunk[k] = unpack_one_tobytes(ct_tobytes, 6*i + k);
			hinv_chunk[k] = unpack_one_tobytes(hinv_tobytes, 6*i + k);
			diff[k] = ct_chunk[k] - a->coeffs[6*i + k];
		}

		basemul(r->coeffs + 6*i,     diff,     hinv_chunk,      zetas[144 + i]);
		basemul(r->coeffs + 6*i + 3, diff + 3, hinv_chunk + 3, -zetas[144 + i]);
	}
}

/* [OPT8] poly_sotp_inv_compressed: sotp_inv from compressed ternary. */
static int16_t decompress_one_ternary(const uint8_t *m1t, int idx)
{
	uint8_t bits = (m1t[idx >> 2] >> (2 * (idx & 3))) & 0x3;
	return (bits == 3) ? -1 : (int16_t)bits;
}

int poly_sotp_inv_compressed(unsigned char *msg, const uint8_t *m1t, const unsigned char *buf)
{
	uint32_t t1, t2, t3, t4;
	uint32_t r = 0;

	for(int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			t1 = load16_littleendian(buf + 32*i + 2*j);
			t2 = load16_littleendian(buf + 32*i + 2*j + 108);
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

	for (int j = 0; j < 4; j++)
	{
		t1 = load16_littleendian(buf + 96 + 2*j);
		t2 = load16_littleendian(buf + 96 + 2*j + 108);
		t3 = 0;

		for(int k = 0; k < 16; k++)
		{
			t4 = t2 & 0x1;
			t4 += decompress_one_ternary(m1t, 768 + 4*k + j);
			r |= t4;
			t4 = (t4^t1) & 0x1;
			t3 ^= t4 << k;

			t1 >>= 1;
			t2 >>= 1;
		}

		msg[96 + 2*j   ] = t3;
		msg[96 + 2*j + 1] = t3 >> 8;
	}

	for (int j = 0; j < 2; j++)
	{
		t1 = load16_littleendian(buf + 104 + 2*j);
		t2 = load16_littleendian(buf + 104 + 2*j + 108);
		t3 = 0;

		for(int k = 0; k < 16; k++)
		{
			t4 = t2 & 0x1;
			t4 += decompress_one_ternary(m1t, 832 + 2*k + j);
			r |= t4;
			t4 = (t4^t1) & 0x1;
			t3 ^= t4 << k;

			t1 >>= 1;
			t2 >>= 1;
		}

		msg[104 + 2*j   ] = t3;
		msg[104 + 2*j + 1] = t3 >> 8;
	}

	r = r >> 1;
	r = (-(uint64_t)r) >> 63;

	return r;
}
