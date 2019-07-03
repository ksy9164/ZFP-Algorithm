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
    for ( int i = 0; i < INPUT_SIZE; i++ ) {
        zfp->original[i] = ((double)(rand()%100000))/1000;
        zfp->decompressed[i] = 0;
        printf( "%04.04lf ", zfp->original[i] );
        if ( (i+1)%4 == 0 ) printf( "\n" );
    }
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

    if (exponent > 0) {
        for (int i = 0; i < exponent; ++i) {
            ans *= base;
        }
        return ans;
    }
    else if (exponent < 0){
        for (int i = 0; i < -exponent; ++i) {
            ans /= base;
        }
        return ans;
    }
    else
        return 1;
}

void convert_to_fixtedpoint(struct zfp_data *zfp)
{

    /** Get maximum exp */
    int exp_max = -EXP_MAX;
    for ( int i = 0; i < INPUT_SIZE; i++ ) {
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
   for ( int i = 0; i < INPUT_SIZE; i++ ) {
       zfp->idata[i] = (int64_t)(zfp->original[i]*(pow(2, 64-2 - exp_max)));
       printf( "%ld ,", zfp->idata[i] );
       if ( (i+1)%4 == 0 ) printf( "\n" );
   }

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
    for ( int i = 0; i < INPUT_SIZE; i++ ) {
        zfp->udata[i] = int2uint_int64(zfp->idata[perm_2[i]]);
        printf( "%016lx ", zfp->udata[i] );
        if ( (i+1)%4 == 0 ) printf( "\n" );
    }
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
         * Ex) INPUT_SIZE element data -> 1 0 0 1 1 0 0 1 1 0 1 0 1 1 1 0 */
        for ( int j = 0; j < INPUT_SIZE; j++ ) {
            bitplane |= (((zfp->udata[j]>>(num_bitplanes-i-1))&1)<<j);
        }

        int verbatim_emit = MIN(verbatim_bits, zfp->bit_budget);

        EncodeBits(zfp, bitplane, verbatim_emit);
        bitplane >>= verbatim_emit;

        total_bits += verbatim_emit;

        zfp->bit_budget -= verbatim_emit;
        if ( zfp->bit_budget <= 0 ) break;

        int left = INPUT_SIZE-verbatim_emit;
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
    printf("---------------------COMPRESSION IS DONE!!---------------------\n");
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

void decompress(struct zfp_data *zfp)
{
    restore_embeded_data(zfp);

    restore_transed_data(zfp);
}

void DecodeBit(struct zfp_data *zfp, uint64_t *dst)
{
    int byteoff = (zfp->decoff/8);
    int bitoff = zfp->decoff%8;
    uint8_t buf = zfp->buffer[byteoff];
    *dst = (buf>>bitoff) & 1;
    zfp->decoff++;
}
void DecodeBits(struct zfp_data *zfp,uint64_t *dst, int bits)
{
    *dst = 0;
    for ( int i = 0; i < bits; i++ ) {
        uint64_t o = 0;
        DecodeBit(zfp, &o);
        *dst |= (o<<i);
    }
}
void restore_embeded_data(struct zfp_data *zfp)
{
    int bit_budget = zfp->curbits;
    uint64_t e;
    DecodeBits(zfp, &e, EBITS);
    zfp->exp_max = ((int)e) - EXP_MAX;
    int precision_max = MIN(ZFP_MAX_PREC, MAX(0, zfp->exp_max - ZFP_MIN_EXP + (2*(zfp->dimension+1))));
    printf("---------------------START DECOMPRESS DATA---------------------\n");
    printf( "e: %ld exp_max retrieved to be %d, precision_max: %d\n", e, zfp->exp_max, precision_max );
 
    zfp->udata_restored = (uint64_t *)malloc(sizeof(uint64_t) * INPUT_SIZE);

    int num_bitplanes = MIN(sizeof(uint64_t)*8, (size_t)precision_max);
    int verbatim_bits = 0;
    int verbatim_emit = 0;
    int bitplane_cnt = 0;

    for ( int i = 0; i < num_bitplanes; i++ ) {
        bitplane_cnt ++;
        int bitplane_emitted_bits = 0;
        int left = INPUT_SIZE;

        uint64_t bitplane = 0;
        int boff = 0;
        verbatim_emit = MIN(verbatim_bits, bit_budget);

        if ( verbatim_bits > 0 ) {
            DecodeBits(zfp, &bitplane, verbatim_emit);
            bit_budget -= verbatim_emit;
            left -= verbatim_emit;
            boff = verbatim_emit;
            bitplane_emitted_bits += verbatim_emit;
        }

        while ( left > 0 && bit_budget > 0 ) {
            uint64_t t;
            DecodeBit(zfp, &t);
            bit_budget --;
            bitplane_emitted_bits++;
            if ( t != 0 ) { // rest is zero
                break;
            } else {
                while ( t == 0 && left > 0 && bit_budget > 0 ) {
                    DecodeBit(zfp, &t);
                    if ( t != 0 ) {
                        bitplane |= ( 1<<boff );
                    }
                    boff++;
                    verbatim_bits++;
                    bitplane_emitted_bits++;
                    bit_budget --;
                    left --;
                }
            }
        }

        for ( int j = 0; j < INPUT_SIZE; j++ ) {
            zfp->udata_restored[j] |= (((bitplane>>j)&1)<<(num_bitplanes-i-1));
        }
        if ( bit_budget <= 0 ) break;
    }

    printf("---------------------RESTORED EMBEDED DATA---------------------\n");
    for ( int i = 0; i < 16; i++ ) {
        printf( "%016lx ", zfp->udata_restored[i] );
        if ( (i+1)%4 == 0 ) printf( "\n" );
    }
}
void restore_transed_data(struct zfp_data *zfp)
{
    int64_t idata[16];
 
    for ( int i = 0; i < 16; i++ )
        idata[perm_2[i]] = uint2int_uint64(zfp->udata_restored[i]);
 
    for (int x = 0; x < 4; x++)
        inv_lift_int64(idata + 1 * x, 4);
 
    for (int y = 0; y < 4; y++)
        inv_lift_int64(idata + 4 * y, 1);
 
    printf("---------------------RESTORED TRANSED DATA---------------------\n");
    for ( int i = 0; i < 16; i++ ) {
        printf( "%ld ,", idata[i] );
        if ( (i+1)%4 == 0 ) printf( "\n" );
    }
    printf("-----------------------DECOMPRESSED DATA-----------------------\n");
    for ( int i = 0; i < 16; i++ ) {
        double q = pow(2,zfp->exp_max-(64-2));
        zfp->decompressed[i] = ((double)idata[i]*q);
        printf( "%04.04lf ", zfp->decompressed[i] );
        if ( (i+1)%4 == 0 ) printf( "\n" );
    }
}
static void inv_lift_int64(int64_t* p, uint64_t s)
{
  int64_t x, y, z, w;
  x = *p; p += s;
  y = *p; p += s;
  z = *p; p += s;
  w = *p; p += s;
 
  y += w >> 1; w -= y >> 1;
  y += w; w <<= 1; w -= y;
  z += x; x <<= 1; x -= z;
  y += z; z <<= 1; z -= y;
  w += x; x <<= 1; x -= w;
 
  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}

void show_loss(struct zfp_data *zfp)
{
    double error = 0;
    printf("---------------------TOTAL COMPRESSION ERR---------------------\n");
    for ( int i = 0; i < 16; i++ ) {
        error += abs_double(zfp->original[i]-zfp->decompressed[i]);
    }
    printf( "Total error: %lf\n", error );
}
double abs_double(double d)
{
    if (d < 0)
        return -d;
    else
        return d;
}

