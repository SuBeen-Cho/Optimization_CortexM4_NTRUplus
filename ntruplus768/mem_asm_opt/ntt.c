#include <stdint.h>
#include "params.h"
#include "ntt.h"

#define NTRUPLUS_R            -147 // R = 2^16 mod q
#define NTRUPLUS_RINV         -682 // (R)^(-1) mod q
#define NTRUPLUS_RSQ           867 // (R^2) mod q
#define NTRUPLUS_QINV        12929 // (q)^(-1) mod (2^16)

#define NTRUPLUS_OMEGA        -886 // (omega * R) mod q
#define NTRUPLUS_ZMINUSZ5INV -1665 // (z - z^5)^(-1) * R mod q
                                   // where z = zeta^((n/d)/6)

#define NTRUPLUS_NINV         -811 // (n/d)^(-1) * R mod q
#define NTRUPLUS_2NINV       -1622 // 2 * (n/d)^(-1) * R mod q

// zetas: Montgomery-form twiddle factors
const int16_t zetas[192] = {
	 -147, -1033,  -682,  -248,  -708,   682,     1,  -722,
	 -723,  -257, -1124,  -867,  -256,  1484,  1262, -1590,
	 1611,   222,  1164, -1346,  1716, -1521,  -357,   395,
	 -455,   639,   502,   655,  -699,   541,    95, -1577,
	-1241,   550,   -44,    39,  -820,  -216,  -121,  -757,
	 -348,   937,   893,   387,  -603,  1713, -1105,  1058,
	 1449,   837,   901,  1637,  -569, -1617, -1530,  1199,
	   50,  -830,  -625,     4,   176,  -156,  1257, -1507,
	 -380,  -606,  1293,   661,  1428, -1580,  -565,  -992,
	  548,  -800,    64,  -371,   961,   641,    87,   630,
	  675,  -834,   205,    54, -1081,  1351,  1413, -1331,
	-1673, -1267, -1558,   281, -1464,  -588,  1015,   436,
	  223,  1138, -1059,  -397,  -183,  1655,   559, -1674,
	  277,   933,  1723,   437, -1514,   242,  1640,   432,
	-1583,   696,   774,  1671,   927,   514,   512,   489,
	  297,   601,  1473,  1130,  1322,   871,   760,  1212,
	 -312,  -352,   443,   943,     8,  1250,  -100,  1660,
	  -31,  1206, -1341, -1247,   444,   235,  1364, -1209,
	  361,   230,   673,   582,  1409,  1501,  1401,   251,
	 1022, -1063,  1053,  1188,   417, -1391,   -27, -1626,
	 1685,  -315,  1408, -1248,   400,   274, -1543,    32,
	-1550,  1531, -1367,  -124,  1458,  1379,  -940, -1681,
	   22,  1709,  -275,  1108,   354, -1728,  -968,   858,
	 1221,  -218,   294,  -732, -1095,   892,  1588,  -779
};

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              a 16-bit integer congruent to a * R^-1 mod q,
*              where R = 2^16.
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           must lie in {-q*2^15, ..., q*2^15-1}
*
* Returns:     an integer in {-q+1, ..., q-1} congruent to
*              a * R^-1 mod q.
**************************************************/
static inline int16_t montgomery_reduce(int32_t a)
{
	int16_t t;
	
	t = (int16_t)a * NTRUPLUS_QINV;
	t = (a - (int32_t)t * NTRUPLUS_Q) >> 16;
	return t;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes a
*              centered representative congruent to a mod q in
*              {-(q+1)/2, ..., (q+1)/2}.
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {-(q+1)/2, ..., (q+1)/2} congruent to a mod q.
**************************************************/
static inline int16_t barrett_reduce(int16_t a)
{
	int16_t t;
	const int16_t v = ((1<<26) + NTRUPLUS_Q/2) / NTRUPLUS_Q;
	
	t  = ((int32_t)v*a + (1<<25)) >> 26;
	t *= NTRUPLUS_Q;
	return a - t;
}

/*************************************************
* Name:        fqmul
*
* Description: Multiplication followed by Montgomery reduction.
*
* Arguments:   - int16_t a: first factor
*              - int16_t b: second factor
*
* Returns:     16-bit integer congruent to a*b*R^-1 mod q.
**************************************************/
static inline int16_t fqmul(int16_t a, int16_t b)
{
    return montgomery_reduce((int32_t)a * b);
}


/* [FIX-ASM] fqinv removed — it was only used by the (now removed) C baseinv
 * fallback. The asm baseinv (ntt_baseinv.S) inlines its own clean-convention
 * fqinv addition chain. */

/*************************************************
* Name:        ntt
*
* Description: Number-theoretic transform (NTT) in R_q. Transforms the
*              coefficient representation of r into a representation
*              where each block of 4 coefficients corresponds to an
*              element of Zq[X]/(X^4 - zeta_i).
*
* Arguments:   - int16_t r[NTRUPLUS_N]: pointer to input/output vector;
*                                       input coefficients of r in R_q,
*                                       output NTT representation in the
*                                       product ring Zq[X]/(X^4 - zeta_i)
*
* Returns:     none.
**************************************************/
static void ntt_inplace(int16_t r[NTRUPLUS_N])   /* [PORT] in-place asm driver */
{
	int16_t t1, t2, t3;
	int16_t zeta1, zeta2;
	int k = 1;

	zeta1 = zetas[k++];

	for (int i = 0; i < NTRUPLUS_N / 2; i++)
	{
		t1 = fqmul(zeta1, r[i + NTRUPLUS_N / 2]);

		r[i + NTRUPLUS_N / 2] = r[i] + r[i + NTRUPLUS_N / 2] - t1;
		r[i                 ] = r[i]                         + t1;
	}

	// radix-3 (OMEGA) layer in asm; consumes 2*(N/384)=4 zetas. KAT-gated.
	ntt_radix3_layer(r, &zetas[k], NTRUPLUS_N, 128);
	k += 4;

	// radix-2 Cooley-Tukey layers (step 64..4) + Barrett pass, done in asm
	// (bit-identical to the C reference; KAT-gated). zetas index k==6 here.
	ntt_ct_layers(r, &zetas[k]);
}

/*************************************************
* Name:        invntt
*
* Description: Inverse number-theoretic transform (NTT) in R_q. Transforms
*              the NTT representation in r, where each block of 4
*              coefficients corresponds to an element of Zq[X]/(X^4 - zeta_i),
*              back to the coefficient representation in R_q.
*
* Arguments:   - int16_t r[NTRUPLUS_N]: pointer to input/output vector;
*                                       input NTT representation in the
*                                       product ring Zq[X]/(X^4 - zeta_i),
*                                       output coefficient representation
*                                       in R_q
*
* Returns:     none.
**************************************************/
static void invntt_inplace(int16_t r[NTRUPLUS_N])   /* [PORT] in-place asm driver */
{
	int16_t t1, t2, t3;
	int16_t zeta1, zeta2;
	int k = 191;

	// GS radix-2 layers (step 4..64) in asm; zeta index k decreases. KAT-gated.
	for (int step = 4; step <= 64; step <<= 1)
	{
		invntt_gs_layer(r, &zetas[k], NTRUPLUS_N, step);
		k -= NTRUPLUS_N / (2 * step);
	}

	// inverse radix-3 (OMEGA) layer in asm (768: 1 layer, no barrett). KAT-gated.
	invntt_radix3_layer(r, &zetas[k], NTRUPLUS_N, 128, 0);
	k -= 4;

	// final X^(n/2) merge + scaling in asm.
	// [FIX] The asm INVNTT layers leave the result one Montgomery factor low
	// (asm = clean * R^-1) vs the clean reference. The asm convention paired
	// this with the (now removed) R^0 asm basemul. To restore clean's domain
	// at the INVNTT output, pass clean's final-scale constants (1679 = ninv*R,
	// -99 = 2ninv*R) instead of (NINV, 2NINV). invntt_scale computes
	// mont(c*(T1-T2)) / mont(c*T2), so this is exactly clean's fqmul(1679,..)/
	// fqmul(-99,..) final loop -> asm INVNTT becomes byte-identical to clean.
	invntt_scale(r, NTRUPLUS_N, 1679, -99);
}

/*************************************************
* [PORT] 2-arg (out,in) wrappers — match the memory-opt poly.c interface.
* The asm drivers are in-place, so we copy a->r first (skip if aliased).
* Host twin verified: ntt/invntt here back the full opt9 memory layer with
* byte-identical KAT (ntt_ternary/invntt_lazy are not needed — see report).
**************************************************/
void ntt(int16_t r[NTRUPLUS_N], const int16_t a[NTRUPLUS_N])
{
	if (r != a)
		for (int i = 0; i < NTRUPLUS_N; i++) r[i] = a[i];
	ntt_inplace(r);
}

void invntt(int16_t r[NTRUPLUS_N], const int16_t a[NTRUPLUS_N])
{
	if (r != a)
		for (int i = 0; i < NTRUPLUS_N; i++) r[i] = a[i];
	invntt_inplace(r);
}

/* [FIX-ASM] baseinv / basemul / basemul_add are now provided by the asm
 * (ntt_baseinv.S / ntt_basemul.S) with the Montgomery convention corrected to
 * clean (basemul/add: removed trailing *RSQ; baseinv: removed trailing *RINV in
 * fqinv). The previous clean-convention C fallbacks are removed so the asm is
 * the single definition. Declarations remain in ntt.h. */
