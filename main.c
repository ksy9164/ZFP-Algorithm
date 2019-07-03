#include "zfp.h"

int main(void) {
    srand(time(NULL));

    struct zfp_data *zfp = (struct zfp_data *)malloc(sizeof(struct zfp_data));
    
    /** Initialize Data */
    init(zfp);

    /** Compress Data */
    compress(zfp);

    /** Decompress Data */
//    decompress(zfp);

    /** Show loss */
//    show_loss(zfp);
}
