#include <cuda_runtime.h>
#include "CUDADefines.h"
#include "tetmesh.h"



__global__ void run_reduction(bool *con, bool *blockCon,int* ActiveList, int nActiveBlock, int* blockSizes)
{
  int list_idx = blockIdx.y*gridDim.x + blockIdx.x;
  int maxblocksize = blockDim.x;
  int tx = threadIdx.x;
  int block_idx = ActiveList[list_idx];

  int blocksize = blockSizes[block_idx];

  __shared__ bool s_block_conv;


  s_block_conv = true;
  __syncthreads();

  if(tx < blocksize)
  {
    if(!con[maxblocksize*block_idx+tx])
      s_block_conv= false;
  }
  __syncthreads();

  if(tx == 0)
  {
    blockCon[block_idx] = s_block_conv; // active list is negation of tile convergence (active = not converged)
  }
}

__device__ bool operator==(const float3 & a, const float3 & b)
{    return (a.x==b.x && a.y==b.y  && a.z==b.z);}

__device__ float3 operator+(const float3 & a, const float3 & b)
{    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);}

__device__ float3 operator-(const float3 & a, const float3 & b)
{    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);}

__device__ float3 operator*(const float3 & a, const float3 & b)
{    return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);}

__device__ float3 operator/(const float3 & a, const float3 & b)
{    return make_float3(a.x/b.x, a.y/b.y, a.z/b.z);}

__device__ float operator ^ (const float3 & a, const float3 & b)
{   return (float)a.x * b.x + a.y*b.y + a.z*b.z;}

__device__ float3 operator % (const float3 & v1, const float3 & v2)
{   return make_float3(v1.y*v2.z - v1.z*v2.y,
    v1.z*v2.x - v1.x*v2.z,
    v1.x*v2.y - v1.y*v2.x);}

__device__ float norm(const float3 & a)
{   return (float)sqrt(a.x*a.x+a.y*a.y+a.z*a.z);}

__device__ void normalize(float3 & a)
{
  float n = norm(a);
  a.x /= n;
  a.y /= n;
  a.z /= n;
}

__device__ float localSolverTet1(float TA, float TB, float TC, float ACAC, float ACBC, float ACCD, float BCBC, float BCCD, float CDCD)
{


  if(TA >= LARGENUM && TB >= LARGENUM && TC >= LARGENUM)
    return LARGENUM;
  float p, q, r;
  float lambda1, lambda2, lambda3;
  float FaceTAC = LARGENUM, FaceTAB = LARGENUM, FaceTBC = LARGENUM;
  float delta, TE;
  float TD = LARGENUM;
  float TAC = TC - TA;
  float TBC = TC - TB;
  float TAB = TB - TA;
  //calculate FaceTBC, let lambda1 = 0
  p = BCBC*TBC*TBC - BCBC*BCBC;
  q = (BCCD+BCCD)*TBC*TBC - 2*BCBC*BCCD;
  r = TBC*TBC*CDCD - BCCD*BCCD;

  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = 0.0;
    lambda2 = (-q + sqrt(delta)) / (2.0*p);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda2 >= 0 && lambda2 <= 1)
    {
      TE = TB*lambda2 + TC*lambda3;
      FaceTBC = min(FaceTBC,TE + sqrt(/*(lambda DOT theta) */(lambda2*BCBC+BCCD)* lambda2 + (lambda2*BCCD+CDCD)));
    }
    lambda1 = 0.0;
    lambda2 = (-q - sqrt(delta)) / (2.0*p);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda2 >= 0 && lambda2 <= 1)
    {
      TE = TB*lambda2  + TC*lambda3;
      FaceTBC = min(FaceTBC, TE + sqrt(/*(lambda DOT theta) */(lambda2*BCBC + BCCD)* lambda2 + (lambda2*BCCD + CDCD)));
    }
  }
  FaceTBC = min(FaceTBC, min(TB + sqrt((BCBC + BCCD) + (BCCD + CDCD)), TC + sqrt(CDCD)));
  float3 gamma = make_float3(ACAC-ACBC, ACBC-BCBC, ACCD-BCCD);
  p = (TAB*TAB*ACAC - gamma.x*gamma.x)  + ( BCBC*TAB*TAB - gamma.y*gamma.y ) - ( (ACBC+ACBC)*TAB*TAB - 2*gamma.x*gamma.y );

  q = -(BCBC*TAB*TAB - gamma.y*gamma.y)*2 +
    ( (ACBC+ACBC)*TAB*TAB - 2*gamma.x*gamma.y ) +
    ( (ACCD+ACCD)*TAB*TAB - 2*gamma.x*gamma.z ) -
    ( (BCCD+BCCD)*TAB*TAB - 2*gamma.y*gamma.z );

  r = (TAB*TAB*BCBC - gamma.y*gamma.y) + ( (BCCD+BCCD)*TAB*TAB - 2*gamma.y*gamma.z ) + (TAB*TAB*CDCD - gamma.z*gamma.z);


  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = (-q + sqrt(delta)) / (2.0*p);
    lambda2 = 1-lambda1;
    lambda3 = 0.0;

    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2;
      FaceTAB = min(FaceTAB, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));


    }

    lambda1 = (-q - sqrt(delta)) / (2.0*p);
    lambda2 = 1-lambda1;
    lambda3 = 0.0;

    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2;
      FaceTAB = min(FaceTAB, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));
    }
  }
  FaceTAB = min(FaceTAB, min(TB + sqrt((BCBC + BCCD) + (BCCD + CDCD)), TA + sqrt((ACAC + ACCD) + (ACCD + CDCD))));
  //calculate FaceTAC, let lambda2 = 0
  p = ACAC*TAC*TAC - ACAC*ACAC;
  q = (ACCD+ACCD)*TAC*TAC - 2*ACAC*ACCD;
  r = TAC*TAC*CDCD - ACCD*ACCD;

  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = (-q + sqrt(delta)) / (2.0*p);
    lambda2 = 0.0;
    lambda3 = 1 - lambda1 - lambda2;
    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TC*lambda3;
      FaceTAC = min(FaceTAC, TE + sqrt((lambda1*ACAC + ACCD) * lambda1 + (lambda1*ACCD + CDCD)));
    }
    lambda1 = (-q - sqrt(delta)) / (2.0*p);
    lambda2 = 0.0;
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2 + TC*lambda3;
      FaceTAC = min(FaceTAC, TE + sqrt((lambda1*ACAC + ACCD) * lambda1 + (lambda1*ACCD + CDCD)));
    }
  }
  FaceTAC = min(FaceTAC, min(TA + sqrt((ACAC + ACCD) + (ACCD + CDCD)), TC + sqrt(CDCD)));

  ////////Done calculating FaceTAC/////////////////////////

  float s = TAC*ACBC - TBC*ACAC;
  float t = TAC*BCCD - TBC*ACCD;
  float h = -(TAC*BCBC - TBC*ACBC);


  p = (TAC*TAC*ACAC- ACAC*ACAC)*h*h  + ( BCBC*TAC*TAC - ACBC*ACBC )*s*s + ( (ACBC+ACBC)*TAC*TAC - 2*ACAC*ACBC )*s*h;

  q = (BCBC*TAC*TAC - ACBC*ACBC)*2*s*t +
    ( (ACBC+ACBC)*TAC*TAC - 2*ACAC*ACBC ) * t*h +
    ( (ACCD+ACCD)*TAC*TAC - 2*ACAC*ACCD) * h*h +
    ( (BCCD+BCCD)*TAC*TAC - 2*ACBC*ACCD ) * s*h;

  r = (TAC*TAC*BCBC - ACBC*ACBC)*t*t + ( (BCCD+BCCD)*TAC*TAC - 2*ACBC*ACCD ) * t*h + (TAC*TAC*CDCD- ACCD*ACCD)*h*h;

  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = (-q + sqrt(delta)) / (2.0*p);
    lambda2 = (s*lambda1 + t) / (h+SMALLNUM);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda1 >= 0 && lambda1 <= 1 && lambda2 >= 0 && lambda2 <= 1 && lambda3 >= 0 && lambda3 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2 + TC*lambda3;
      TD = min(TD, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));
    }

    lambda1 = (-q - sqrt(delta)) / (2.0*p);
    lambda2 = (s*lambda1 + t) / (h+SMALLNUM);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda1 >= 0 && lambda1 <= 1 && lambda2 >= 0 && lambda2 <= 1 && lambda3 >= 0 && lambda3 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2 + TC*lambda3;
      TD = min(TD, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));
    }
    TD = min(TD, min(FaceTBC, min(FaceTAB, FaceTAC)));
  }
  else
  {
    TD = min(TD, min(FaceTBC, min(FaceTAB, FaceTAC)));
  }
  return TD;
}

__device__ float localSolverTet2(float TA, float TB, float TC, float ACAC, float ACBC, float ACCD, float BCBC, float BCCD, float CDCD)
{
  if(TA >= LARGENUM && TB >= LARGENUM && TC >= LARGENUM)
    return LARGENUM;
  float p, q, r;
  float lambda1, lambda2, lambda3;
  float FaceTAC = LARGENUM, FaceTAB = LARGENUM, FaceTBC = LARGENUM;
  float delta, TE;
  float TD = LARGENUM;
  float TAC = TC - TA;
  float TBC = TC - TB;
  float TAB = TB - TA;
  //calculate FaceTBC, let lambda1 = 0
  p = BCBC*TBC*TBC - BCBC*BCBC;
  q = (BCCD+BCCD)*TBC*TBC - 2*BCBC*BCCD;
  r = TBC*TBC*CDCD - BCCD*BCCD;

  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = 0.0;
    lambda2 = (-q + sqrt(delta)) / (2.0*p);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda2 >= 0 && lambda2 <= 1)
    {
      TE = TB*lambda2 + TC*lambda3;
      FaceTBC = min(FaceTBC, TE + sqrt(/*(lambda DOT theta) */(lambda2*BCBC + BCCD)* lambda2 + (lambda2*BCCD + CDCD)));
    }
    lambda1 = 0.0;
    lambda2 = (-q - sqrt(delta)) / (2.0*p);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda2 >= 0 && lambda2 <= 1)
    {
      TE = TB*lambda2  + TC*lambda3;
      FaceTBC = min(FaceTBC, TE + sqrt(/*(lambda DOT theta) */(lambda2*BCBC + BCCD)* lambda2 + (lambda2*BCCD + CDCD)));
    }
  }
  FaceTBC = min(FaceTBC, min(TB + sqrt((BCBC + BCCD) + (BCCD + CDCD)), TC + sqrt(CDCD)));
  //calculate FaceTAB, let lambda3 = 0

  float3 gamma = make_float3(ACAC-ACBC, ACBC-BCBC, ACCD-BCCD);
  p = (TAB*TAB*ACAC - gamma.x*gamma.x)  + ( BCBC*TAB*TAB - gamma.y*gamma.y ) - ( (ACBC+ACBC)*TAB*TAB - 2*gamma.x*gamma.y );

  q = -(BCBC*TAB*TAB - gamma.y*gamma.y)*2 +
    ( (ACBC+ACBC)*TAB*TAB - 2*gamma.x*gamma.y ) +
    ( (ACCD+ACCD)*TAB*TAB - 2*gamma.x*gamma.z ) -
    ( (BCCD+BCCD)*TAB*TAB - 2*gamma.y*gamma.z );

  r = (TAB*TAB*BCBC - gamma.y*gamma.y) + ( (BCCD+BCCD)*TAB*TAB - 2*gamma.y*gamma.z ) + (TAB*TAB*CDCD - gamma.z*gamma.z);


  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = (-q + sqrt(delta)) / (2.0*p);
    lambda2 = 1-lambda1;
    lambda3 = 0.0;
    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2;
      FaceTAB = min(FaceTAB, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));
    }
    lambda1 = (-q - sqrt(delta)) / (2.0*p);
    lambda2 = 1-lambda1;
    lambda3 = 0.0;

    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2;
      FaceTAB = min(FaceTAB, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));
    }
  }
  FaceTAB = min(FaceTAB, min(TB + sqrt((BCBC + BCCD) + (BCCD + CDCD)), TA + sqrt((ACAC + ACCD) + (ACCD + CDCD))));
  //calculate FaceTAC, let lambda2 = 0
  p = ACAC*TAC*TAC - ACAC*ACAC;
  q = (ACCD+ACCD)*TAC*TAC - 2*ACAC*ACCD;
  r = TAC*TAC*CDCD - ACCD*ACCD;

  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = (-q + sqrt(delta)) / (2.0*p);
    lambda2 = 0.0;
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TC*lambda3;
      FaceTAC = min(FaceTAC, TE + sqrt((lambda1*ACAC + ACCD) * lambda1 + (lambda1*ACCD + CDCD)));
    }

    lambda1 = (-q - sqrt(delta)) / (2.0*p);
    lambda2 = 0.0;
    lambda3 = 1 - lambda1 - lambda2;
    if(lambda1 >= 0 && lambda1 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2 + TC*lambda3;
      FaceTAC = min(FaceTAC, TE + sqrt((lambda1*ACAC + ACCD) * lambda1 + (lambda1*ACCD + CDCD)));
    }
  }
  FaceTAC = min(FaceTAC, min(TA + sqrt((ACAC + ACCD) + (ACCD + CDCD)), TC + sqrt(CDCD)));

  ////////Done calculating FaceTAC/////////////////////////

  float s = TAC*ACBC - TBC*ACAC;
  float t = TAC*BCCD - TBC*ACCD;
  float h = -(TAC*BCBC - TBC*ACBC);


  p = (TAC*TAC*ACAC- ACAC*ACAC)*h*h  + ( BCBC*TAC*TAC - ACBC*ACBC )*s*s + ( (ACBC+ACBC)*TAC*TAC - 2*ACAC*ACBC )*s*h;

  q = (BCBC*TAC*TAC - ACBC*ACBC)*2*s*t +
    ( (ACBC+ACBC)*TAC*TAC - 2*ACAC*ACBC ) * t*h +
    ( (ACCD+ACCD)*TAC*TAC - 2*ACAC*ACCD) * h*h +
    ( (BCCD+BCCD)*TAC*TAC - 2*ACBC*ACCD ) * s*h;

  r = (TAC*TAC*BCBC - ACBC*ACBC)*t*t + ( (BCCD+BCCD)*TAC*TAC - 2*ACBC*ACCD ) * t*h + (TAC*TAC*CDCD- ACCD*ACCD)*h*h;

  delta = q*q - 4*p*r;

  if(delta >= 0.0)
  {
    lambda1 = (-q + sqrt(delta)) / (2.0*p);
    lambda2 = (s*lambda1 + t) / (h+SMALLNUM);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda1 >= 0 && lambda1 <= 1 && lambda2 >= 0 && lambda2 <= 1 && lambda3 >= 0 && lambda3 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2 + TC*lambda3;
      TD = min(TD, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));
    }
    lambda1 = (-q - sqrt(delta)) / (2.0*p);
    lambda2 = (s*lambda1 + t) / (h+SMALLNUM);
    lambda3 = 1 - lambda1 - lambda2;

    if(lambda1 >= 0 && lambda1 <= 1 && lambda2 >= 0 && lambda2 <= 1 && lambda3 >= 0 && lambda3 <= 1)
    {
      TE = TA*lambda1 + TB*lambda2 + TC*lambda3;
      TD = min(TD, TE + sqrt((lambda1*ACAC + lambda2*ACBC + ACCD) * lambda1 + (lambda1*ACBC + lambda2*BCBC + BCCD) * lambda2 + (lambda1*ACCD + lambda2*BCCD + CDCD)));
    }
    TD = min(TD, min(FaceTBC, min(FaceTAB, FaceTAC)));
  }
  else
  {
    TD = min(TD, min(FaceTBC, min(FaceTAB, FaceTAC)));
  }
  return TD;
}
extern __shared__ char s_array[];


__global__ void FIMCuda(float3* d_tetMem0,float3* d_tetMem1, float4* d_tetT, float* d_vertT, float* d_speedInv, int* d_vertMem, int* d_vertMemOutside, int* d_BlockSizes, bool* d_con, int* d_ActiveList, int m_maxNumInVert, int m_maxVertMappingInside, int m_maxVertMappingOutside, int nIter)
{
  int list_idx = blockIdx.y*gridDim.x + blockIdx.x;
  // retrieve actual block index from the active list
  int block_idx = d_ActiveList[list_idx];
  int block_size = d_BlockSizes[block_idx];
  int maxNumTotalTets = blockDim.x;
  ///////////////////initialize shared memory//////////////////////////////////////////
  int tx = threadIdx.x;
  int tet_base = block_idx*maxNumTotalTets;
  int vert_base = block_idx*m_maxNumInVert;

  float*      s_tetT = (float*)s_array;
  int*   s_vertMem = (int*)&s_tetT[maxNumTotalTets * 4];

  float4 tetT = d_tetT[tet_base + tx];

  s_tetT[tx*4 + 0] = tetT.x;
  s_tetT[tx*4 + 1] = tetT.y;
  s_tetT[tx*4 + 2] = tetT.z;
  s_tetT[tx*4 + 3] = tetT.w;

  if(tx < m_maxNumInVert)
    for(int i = 0; i< m_maxVertMappingInside; i++)
    {
      s_vertMem[tx*m_maxVertMappingInside + i] = d_vertMem[vert_base*m_maxVertMappingInside + m_maxVertMappingInside*tx + i] - tet_base*4;

    }

  __syncthreads();

  ////////////////done shared memory copy//////////////////////////////////////////////////

  float oldT, newT;
  float TD, TA, TB, TC;

  float3 tetmem0 = d_tetMem0[tet_base + tx];
  float3 tetmem1 = d_tetMem1[tet_base + tx];

  float ACBC = tetmem0.x;
  float BCCD = tetmem0.y;
  float ACCD = tetmem0.z;
  float ADBD = tetmem1.x;
  float ACAD = tetmem1.y;
  float BCBD = tetmem1.z;

  float ADBC = ACBC + BCCD;
  float ACBD = ACBC + ACCD;
  float CDBD = ADBD - ACBD;
  float ADCD = ADBD - ADBC;
  float CDCD = ADCD - ACCD;
  float ADAD = ADCD + ACAD;
  float ACAC = ACAD - ACCD;
  float BDBD = BCBD + CDBD;
  float BCBC = BCBD - BCCD;

  for(int iter=0; iter< nIter; iter++)
  {
    if(tx < block_size)
    {
      oldT = s_tetT[s_vertMem[tx*m_maxVertMappingInside+0]];
      TD = oldT;
    }
    __syncthreads();

    TA = s_tetT[tx*4+0];
    TB = s_tetT[tx*4+1];
    TC = s_tetT[tx*4+2];
    TD = s_tetT[tx*4+3];

    s_tetT[tx * 4 + 3] = min(TD, localSolverTet1(TA, TB, TC, ACAC, ACBC, ACCD, BCBC, BCCD, CDCD));
    s_tetT[tx * 4 + 0] = min(TA, localSolverTet1(TB, TD, TC, BCBC, -BCCD, -ACBC, CDCD, ACCD, ACAC));
    s_tetT[tx * 4 + 1] = min(TB, localSolverTet1(TA, TD, TC, ACAC, -ACCD, -ACBC, CDCD, BCCD, BCBC));
    s_tetT[tx * 4 + 2] = min(TC, localSolverTet1(TA, TB, TD, ADAD, ADBD, -ADCD, BDBD, -CDBD, CDCD));

    __syncthreads();

    newT = LARGENUM;

    if(tx < block_size)    //block_size is the vertices in this block and it is about warp size so there is no severe divergence
    {
      int tmp2;
      for(int j = 0; (j < m_maxVertMappingInside) && ( (tmp2 = s_vertMem[tx * m_maxVertMappingInside + j]) > -1); j++)                        // find the min and keep the old T for convergence check
      {
        newT = min(newT, s_tetT[tmp2]);
      }

      for(int j =0; (j < m_maxVertMappingInside) && ( (tmp2 = s_vertMem[tx * m_maxVertMappingInside + j]) > -1 ); j++) // update all the old to the min
      {
        s_tetT[tmp2] = newT;
      }
    }
    __syncthreads();
  }

  __syncthreads();

  if(tx < block_size)
  {
    float residue = oldT - newT;
    if (residue < 0.) residue *= -1.;
    d_con[vert_base + tx] = (residue < EPS) ? true : false;
    d_vertT[vert_base + tx] = newT;
  }
  __syncthreads();

}


extern __shared__ float s_run_check_neghbor_array[];
__global__ void run_check_neghbor(float3* d_tetMem0,float3* d_tetMem1, float4* d_tetT, float* d_speedInv, int* d_vertMem, int* d_vertMemOutside,
    int* d_BlockSizes, bool* d_con, int* d_ActiveList,
    int m_maxNumInVert, int m_maxVertMappingInside, int m_maxNumOutVertMapping)
{
  int list_idx = blockIdx.y*gridDim.x + blockIdx.x;

  // retrieve actual block index from the active list
  int block_idx = d_ActiveList[list_idx];
  int block_size = d_BlockSizes[block_idx];
  int maxNumTotalTets = blockDim.x;
  //////////////////initialize shared memory//////////////////////////////////////////

  int tx = threadIdx.x;
  int tet_base = block_idx*maxNumTotalTets;
  int vert_base = block_idx*m_maxNumInVert;

  float*      s_tetT = (float*)s_array;
  int*   s_vertMem = (int*)&s_tetT[maxNumTotalTets * 4];

  float4 tetT = d_tetT[tet_base + tx];

  s_tetT[tx*4 + 0] = tetT.x;
  s_tetT[tx*4 + 1] = tetT.y;
  s_tetT[tx*4 + 2] = tetT.z;
  s_tetT[tx*4 + 3] = tetT.w;

  if(tx < m_maxNumInVert)
    for(int i = 0; i< m_maxVertMappingInside; i++)
    {
      s_vertMem[tx*m_maxVertMappingInside + i] = d_vertMem[vert_base*m_maxVertMappingInside + m_maxVertMappingInside*tx + i] - tet_base*4;

    }

  __syncthreads();

  /////////////done shared memory copy//////////////////////////////////////////////////

  float oldT, newT;
  float TD, TA, TB, TC;

  float3 tetmem0 = d_tetMem0[tet_base + tx];
  float3 tetmem1 = d_tetMem1[tet_base + tx];

  float ACBC = tetmem0.x;
  float BCCD = tetmem0.y;
  float ACCD = tetmem0.z;
  float ADBD = tetmem1.x;
  float ACAD = tetmem1.y;
  float BCBD = tetmem1.z;

  float ADBC = ACBC + BCCD;
  float ACBD = ACBC + ACCD;
  float CDBD = ADBD - ACBD;
  float ADCD = ADBD - ADBC;
  float CDCD = ADCD - ACCD;
  float ADAD = ADCD + ACAD;
  float ACAC = ACAD - ACCD;
  float BDBD = BCBD + CDBD;
  float BCBC = BCBD - BCCD;
  for(int iter=0; iter< 1; iter++)
  {
    if(tx < block_size)
    {
      oldT = s_tetT[s_vertMem[tx*m_maxVertMappingInside+0]];
      TD = oldT;
    }

    __syncthreads();

    TA = s_tetT[tx*4+0];
    TB = s_tetT[tx*4+1];
    TC = s_tetT[tx*4+2];
    TD = s_tetT[tx*4+3];

    s_tetT[tx * 4 + 3] = min(TD, localSolverTet1(TA, TB, TC, ACAC, ACBC, ACCD, BCBC, BCCD, CDCD));
    s_tetT[tx * 4 + 0] = min(TA, localSolverTet1(TB, TD, TC, BCBC, -BCCD, -ACBC, CDCD, ACCD, ACAC));
    s_tetT[tx * 4 + 1] = min(TB, localSolverTet1(TA, TD, TC, ACAC, -ACCD, -ACBC, CDCD, BCCD, BCBC));
    s_tetT[tx * 4 + 2] = min(TC, localSolverTet1(TA, TB, TD, ADAD, ADBD, -ADCD, BDBD, -CDBD, CDCD));

    __syncthreads();

    newT = LARGENUM;

    if(tx < block_size)    //block_size is the vertices in this block and it is about warp size so there is no severe divergence
    {
      int tmp2;
      for(int j = 0; (j < m_maxVertMappingInside) && ( (tmp2 = s_vertMem[tx * m_maxVertMappingInside + j]) > -1); j++)                        // find the min and keep the old T for convergence check
      {
        newT = min(newT, s_tetT[tmp2]);
      }

      for(int j =0; (j < m_maxVertMappingInside) && ( (tmp2 = s_vertMem[tx * m_maxVertMappingInside + j]) > -1 ); j++) // update all the old to the min
      {
        s_tetT[tmp2] = newT;
      }
    }

    __syncthreads();
  }

  __syncthreads();

  if(tx < block_size)
  {
    float residue = oldT - newT;
    d_con[vert_base + tx] = (residue < EPS) ? true : false;
  }
  __syncthreads();

}

__global__ void CopyOutBack(float4* d_tetT, float* d_vertT, int* d_vertMem, int* d_vertMemOutside, int* d_BlockSizes, int* d_ActiveList,int m_maxNumInVert,int m_maxNumTotalTets, int m_maxVertMappingInside, int m_maxVertMappingOutside)
{
  int list_idx = blockIdx.y*gridDim.x + blockIdx.x;
  // retrieve actual block index from the active list
  int block_idx = d_ActiveList[list_idx];
  int block_size = d_BlockSizes[block_idx];

  ////////////initialize shared memory//////////////////////////////////////////

  int tx = threadIdx.x;
  int tet_base = block_idx*m_maxNumTotalTets;
  int vert_base = block_idx*m_maxNumInVert;
  int tmpindex;

  if(tx < block_size)
  {
    float T = d_vertT[vert_base + tx];

    int j =0;
    tmpindex = d_vertMem[block_idx*m_maxVertMappingInside*m_maxNumInVert +  tx * m_maxVertMappingInside + j];
    while(j < m_maxVertMappingInside && (tmpindex - tet_base*4) > -1) // update gloal memory inside all the old to the min
    {
      int segment = tmpindex / 4;
      int offset = tmpindex % 4;
      switch(offset)
      {
      case 0:
        d_tetT[segment].x = T;
        break;
      case 1:
        d_tetT[segment].y = T;
        break;
      case 2:
        d_tetT[segment].z = T;
        break;
      case 3:
        d_tetT[segment].w = T;
        break;
      }
      j++;
      tmpindex = d_vertMem[block_idx*m_maxVertMappingInside*m_maxNumInVert +  tx * m_maxVertMappingInside + j];
    }

    for(int j = 0; (j < m_maxVertMappingOutside) && ((tmpindex = d_vertMemOutside[block_idx*m_maxVertMappingOutside*m_maxNumInVert +  tx * m_maxVertMappingOutside + j]) > -1 ); j++) // update gloal memory outside all the old to the min
    {
      int segment = tmpindex / 4;
      int offset = tmpindex % 4;
      switch(offset)
      {
      case 0:
        d_tetT[segment].x = T;
        break;
      case 1:
        d_tetT[segment].y = T;
        break;
      case 2:
        d_tetT[segment].z = T;
        break;
      case 3:
        d_tetT[segment].w = T;
        break;
      }
    }
  }
}
