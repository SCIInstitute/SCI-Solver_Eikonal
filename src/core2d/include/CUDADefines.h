#ifndef __CUDA_DEFINES_H__
#define __CUDA_DEFINES_H__
#define WARPSIZE 32
#define DIMENSION      3

#define MARK      2
#define SEED      3
#define ACTIVE    1
#define FARP      0
#define STOP      4
#define ALIVE     5
#define TOBEALIVE 6
#define TRIMEMLENGTH   3
#define EDGEMEMLENGTH  3
#define VERTMEMLENGTH  6    //must be even number and best be odd*2 to avoid shared memory bank conflict; 14 for dragon.ts, 12 for dragon.ts_maxSF0.5
#define VERTMEMLENGTHOUTSIDE 6
#define REDUCTIONSHARESIZE 64  //must be power of 2
#define NITER 7   // best 7 for dragon.ts and 7 for dragon.ts_maxSF0.5 one,7 for dragon.ts_maxSF0.5 curvature,


#define EPS (float)1e-06

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


typedef unsigned int uint;
typedef unsigned char uchar;
#endif
