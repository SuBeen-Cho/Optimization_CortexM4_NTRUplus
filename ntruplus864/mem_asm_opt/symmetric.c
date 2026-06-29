// SPDX-License-Identifier: MIT

#include <string.h>
#include "sha2.h"
#include "symmetric.h"
#include "fips202.h"
#include "params.h"

/*
 * [OPT2] All hash functions converted to incremental API.
 * Eliminates large data[] temporary buffers on the stack.
 */

#define HASH_F_MID_BLOCKS ((NTRUPLUS_POLYBYTES - 63) / 64)
#define HASH_F_REMAINDER  ((NTRUPLUS_POLYBYTES - 63) % 64)

void hash_f(uint8_t *buf, const uint8_t *msg)
{
	sha256ctx ctx;
	uint8_t first_block[64];

	sha256_inc_init(&ctx);

	first_block[0] = 0x0;
	memcpy(first_block + 1, msg, 63);
	sha256_inc_blocks(&ctx, first_block, 1);

	sha256_inc_blocks(&ctx, msg + 63, HASH_F_MID_BLOCKS);
	sha256_inc_finalize(buf, &ctx, msg + 63 + HASH_F_MID_BLOCKS * 64, HASH_F_REMAINDER);
}

void hash_g(uint8_t *buf, const uint8_t *msg)
{
	shake256incctx ctx;
	uint8_t domain = 0x1;

	shake256_inc_init(&ctx);
	shake256_inc_absorb(&ctx, &domain, 1);
	shake256_inc_absorb(&ctx, msg, NTRUPLUS_POLYBYTES);
	shake256_inc_finalize(&ctx);
	shake256_inc_squeeze(buf, NTRUPLUS_N / 4, &ctx);
	shake256_inc_ctx_release(&ctx);
}

void hash_h_kem(uint8_t *buf, const uint8_t *msg)
{
	shake256incctx ctx;
	uint8_t domain = 0x2;

	shake256_inc_init(&ctx);
	shake256_inc_absorb(&ctx, &domain, 1);
	shake256_inc_absorb(&ctx, msg, NTRUPLUS_N / 8 + NTRUPLUS_SYMBYTES);
	shake256_inc_finalize(&ctx);
	shake256_inc_squeeze(buf, NTRUPLUS_SSBYTES + NTRUPLUS_N / 4, &ctx);
	shake256_inc_ctx_release(&ctx);
}
