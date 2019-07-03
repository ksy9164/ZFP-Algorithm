#include "zfp.h"

void init(struct zfp_data *zfp)
{
    int bytes = (int)(sizeof(double) * INPUT_SIZE);
    zfp->dimension = DIMENSION;
    zfp->bit_budget = BIT_BUDGET;

    /* Original Input Data & ZFP Decompressed Data */
    zfp->original = (double *)malloc(bytes);
    zfp->decompressed = (double *)malloc(bytes);

    /* Bytes buffer */
    zfp->buffer = (uint8_t *)malloc(bytes);

    for ( int i = 0; i < bytes; i++ )
        zfp->buffer[i] = 0;
    
    zfp->bufferbytes = bytes;
    zfp->curbits = 0;
    zfp->decoff = 0;

    /* Make randomize & init decompressed data */
    printf("------------------------ORIGINAL DATA------------------------\n");
    for ( int i = 0; i < 16; i++ ) {
        zfp->original[i] = ((double)(rand()%100000))/1000;
        zfp->decompressed[i] = 0;
        printf( "%04.04lf ", zfp->original[i] );
        if ( (i+1)%4 == 0 ) printf( "\n" );
    }
    printf("-------------------------------------------------------------\n");
}

int get_exp(double d)
{
    int exp = 0;
    while (d >= 1.0) {
        exp++;
        d /= 2;
    }
    return exp;
}

double pow (double base, double exponent)
{
    double ans = 1;

    for (int i = 0; i < exponent; ++i) {
        ans *= base;
    }

    return ans;
}

void convert_to_fixtedpoint(struct zfp_data *zfp)
{

    /** Get maximum exp */
    int exp_max = -EXP_MAX;
    for ( int i = 0; i < 16; i++ ) {
        if ( zfp->original[i] != 0 ) {
            int exp = get_exp(zfp->original[i]);
            if ( exp > exp_max ) exp_max = exp;
        }
    }

    zfp->precision_max = MIN(ZFP_MAX_PREC, MAX(0, exp_max - ZFP_MIN_EXP + (2*(zfp->dimension+1)))); // (7 + 1074 + 2 * (2 +1) )
    /** If there is no precision */
   if (zfp->precision_max == 0)
       return;

    /** Block biggest exp + EXP_MAX */
   zfp->e = exp_max + EXP_MAX;

   /** write ebits of e */
   printf( "exp_max: %d e: %d\n", exp_max, zfp->e );

   /** Make fixed point value */
   zfp->idata = (int64_t *)malloc(sizeof(int64_t) * INPUT_SIZE);
   printf("------------------------FIXED POINT DATA----------------------\n");
   for ( int i = 0; i < 16; i++ ) {
       zfp->idata[i] = (int64_t)(zfp->original[i]*(pow(2, 64-2 - exp_max)));
       printf( "%ld ,", zfp->idata[i] );
       if ( (i+1)%4 == 0 ) printf( "\n" );
   }
   printf("--------------------------------------------------------------\n");

   return;
}

void block_trans(struct zfp_data *zfp)
{
   if (zfp->precision_max == 0)
       return;

    /** Lifting */
    for (int y = 0; y < 4; y++)
        fwd_lift_int64(zfp->idata + 4 * y, 1);
    for (int x = 0; x < 4; x++)
        fwd_lift_int64(zfp->idata + 1 * x, 4);

    /** Permutate and transfer int to unsigned int */
    zfp->udata = (uint64_t *)malloc(sizeof(uint64_t) * INPUT_SIZE);
   printf("----------------------PERMUTATED UINT DATA---------------------\n");
    for ( int i = 0; i < 16; i++ ) {
        zfp->udata[i] = int2uint_int64(zfp->idata[perm_2[i]]);
        printf( "%016lx ", zfp->udata[i] );
        if ( (i+1)%4 == 0 ) printf( "\n" );
    }
   printf("---------------------------------------------------------------\n");
}

static void fwd_lift_int64(int64_t* p, unsigned int s)
{
  int64_t x, y, z, w;
  x = *p; p += s;
  y = *p; p += s;
  z = *p; p += s;
  w = *p; p += s;
 
  x += w; x >>= 1; w -= x;
  z += y; z >>= 1; y -= z;
  x += z; x >>= 1; z -= x;
  w += y; w >>= 1; y -= w;
  w += y >> 1; y -= w >> 1;
 
  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}

// 0xaaaaaaaaaaaaaaaaUL => 101010101010.... 
// convert negabinary uint to int
static int64_t uint2int_uint64(uint64_t x)
{
  return (int64_t)((x ^ 0xaaaaaaaaaaaaaaaaUL) - 0xaaaaaaaaaaaaaaaaUL);
}

// convert int to negabinary uint // xor (to get negative value)
static uint64_t int2uint_int64(int64_t x)
{
  return ((uint64_t)x + 0xaaaaaaaaaaaaaaaaUL) ^ 0xaaaaaaaaaaaaaaaaUL;
}

void EncodeBit(struct zfp_data *zfp, uint64_t src) {
    int byteoff = (zfp->curbits/8);
    int bitoff = zfp->curbits%8;
    zfp->buffer[byteoff] |= ((src&1)<<bitoff);
 
    zfp->curbits++;
}
void EncodeBits(struct zfp_data *zfp, uint64_t src, int bits) {
    for ( int i = 0; i < bits; i++ ) {
        EncodeBit(zfp, src);
        src >>= 1;
    }
}

void embedded_coding(struct zfp_data *zfp)
{ 

    int total_bits = EBITS; // EBITS = 11 + 1
    if(zfp->precision_max == 0) { // everything is zero...
        printf( "Everything is zero\n" );
        EncodeBit(zfp, 0);
        zfp->bit_budget --;
        total_bits ++;
        return;
    }

    zfp->bit_budget -= EBITS;
    
    EncodeBits(zfp, zfp->e, EBITS);

    // how many lsbs are non-encoded
    int verbatim_bits = 0;

    // for each bit plane
    int bitplane_cnt = 0;
    int num_bitplanes = MIN(sizeof(uint64_t)*8, (size_t)(zfp->precision_max)); // 64

    for ( int i = 0; i < num_bitplanes; i++ ) {
        bitplane_cnt ++;

        uint64_t bitplane = 0;
        
        /** Make INPUTSIZE bit data from each Matrix's elements
         * Get each element's one bit (index is i) data
         * Ex) 16 element data -> 1 0 0 1 1 0 0 1 1 0 1 0 1 1 1 0 */
        for ( int j = 0; j < 16; j++ ) {
            bitplane |= (((zfp->udata[j]>>(num_bitplanes-i-1))&1)<<j);
        }

        int verbatim_emit = MIN(verbatim_bits, zfp->bit_budget);

        EncodeBits(zfp, bitplane, verbatim_emit);
        bitplane >>= verbatim_emit;

        total_bits += verbatim_emit;

        zfp->bit_budget -= verbatim_emit;
        if ( zfp->bit_budget <= 0 ) break;

        int left = 16-verbatim_emit;
        while( left > 0 && zfp->bit_budget > 0 ) {
            if ( bitplane == 0 ) {
                EncodeBit(zfp, 1);

                zfp->bit_budget --;

                total_bits ++;
                break;
            } else {
                EncodeBit(zfp, 0);

                zfp->bit_budget --;
                //verbatim_bits ++;

                total_bits ++;
                while ( (bitplane & 1) != 1 && zfp->bit_budget > 0 && left > 0 ) {
                    EncodeBit(zfp, 0);
                    bitplane >>= 1;

                    left --;
                    zfp->bit_budget --;
                    verbatim_bits ++;

                    total_bits ++;
                }
                if ( (bitplane & 1) == 1 && zfp->bit_budget > 0 && left > 0 ) {
                    EncodeBit(zfp, 1);
                    bitplane >>= 1;

                    left --;
                    zfp->bit_budget --;
                    verbatim_bits ++;

                    total_bits ++;
                }
            }
        }
    }
    printf( "Compression done! emitted %d bits from %ld bits(original) across %d bitplanes\n", total_bits, sizeof(double)*8*INPUT_SIZE, bitplane_cnt );


    printf("Compressed Data is \" ");
    for (int i = 0; i < total_bits / 8; ++i) {
        printf("%" PRIu8 " " ,zfp->buffer[i]);
    }
    printf("\"\n");

}

void compress(struct zfp_data *zfp)
{
    /** Compression algorithm is divided into three steps.
     * 1. Conversion to Fixed-Point
     * 2. Block Transform
     * 3. Embedded Coding */

    convert_to_fixtedpoint(zfp);
    
    block_trans(zfp);

    embedded_coding(zfp);

}
