#ifndef _H_ZFP
#define _H_ZFP

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/* exponent of double is 11 bit signed integer */
#define EXP_MAX ((1<<(11-1))-1)
 
/* exponent of the minimum granularity value expressible via double */
#define ZFP_MIN_EXP -1074
#define ZFP_MAX_PREC 64
#define EBITS (11+1)

/* Size of Matrix */
#define INPUT_SIZE 16
#define DIMENSION 2

/* Set bit budget */
#define BIT_BUDGET 8*8*2

struct zfp_data {
    uint8_t* buffer;
    int bufferbytes;
    int curbits;
    int decoff;

    double *original;
    double *decompressed;

    int dimension;
    int bit_budget;

    int precision_max;

    int e;

    int64_t *idata;
    uint64_t *udata;
};

/* Permutation array */
static const uint8_t perm_2[16] __attribute__((aligned(0x100))) = {
  ((0) + 4 * (0)), // 0  

  ((1) + 4 * (0)), // 1
  ((0) + 4 * (1)), // 4
 
  ((1) + 4 * (1)), // 5
 
  ((2) + 4 * (0)), // 2
  ((0) + 4 * (2)), // 8
 
  ((2) + 4 * (1)), // 6
  ((1) + 4 * (2)), // 9
 
  ((3) + 4 * (0)), // 3
  ((0) + 4 * (3)), // 12
 
  ((2) + 4 * (2)), // 10
 
  ((3) + 4 * (1)), // 7
  ((1) + 4 * (3)), // 13
 
  ((3) + 4 * (2)), // 11
  ((2) + 4 * (3)), // 14
 
  ((3) + 4 * (3)), // 15
};

void init(struct zfp_data *zfp);

void compress(struct zfp_data *zfp);
void convert_to_fixtedpoint(struct zfp_data *zfp);
int get_exp(double d);
double pow (double base, double exponent);
void block_trans(struct zfp_data *zfp);
void embedded_coding(struct zfp_data *zfp);

/* void decompress(struct zfp_data *zfp);
 * void show_loss(struct zfp_data *zfp); */

/* Lift Algorithm */
static void fwd_lift_int64(int64_t* p, unsigned int s);

static int64_t uint2int_uint64(uint64_t x);
static uint64_t int2uint_int64(int64_t x);

void EncodeBit(struct zfp_data *zfp, uint64_t src);
void EncodeBits(struct zfp_data *zfp, uint64_t src, int bits);

#endif
