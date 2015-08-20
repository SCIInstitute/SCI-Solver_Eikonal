
#define WARPSIZE 32
#define DIMENSION      3

#define MARK      2
#define SEED      3
#define ACTIVE    1
#define FARP      0
#define STOP      4
#define ALIVE     5
#define TOBEALIVE 6
#define TETMEMLENGTH   4
#define EDGEMEMLENGTH  3
#define VERTMEMLENGTH  6    //must be even number and best be odd*2 to avoid shared memory bank conflict; 14 for dragon.ts, 12 for dragon.ts_maxSF0.5
#define VERTMEMLENGTHOUTSIDE 6
#define REDUCTIONSHARESIZE 16  //must be power of 2
#define NITER 2   // best 7 for dragon.ts and 7 for dragon.ts_maxSF0.5 one,7 for dragon.ts_maxSF0.5 curvature,


#define EPS 1e-5

typedef unsigned int uint;
typedef unsigned char uchar;
