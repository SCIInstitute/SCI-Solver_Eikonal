#include <cuda_runtime.h>
#include <cutil.h>
#include "CUDADefines.h"
#include "TriMesh.h"



__global__ void run_reduction(int *con, int *blockCon,int* ActiveList, int nActiveBlock, int* blockSizes)
{
  int list_idx = blockIdx.x;
  int tx = threadIdx.x;
  int block_idx = ActiveList[list_idx];
  int start = block_idx*blockDim.x * 2;
  int blocksize = blockSizes[block_idx];
  __shared__ int s_block_conv;
  s_block_conv = 1;
  __syncthreads();

  if (tx < blocksize)
  {
    if (!con[start + tx])
      s_block_conv = 0;
  }
  __syncthreads();

  if(tx == 0)
  {
    blockCon[block_idx] = s_block_conv; // active list is negation of tile convergence (active = not converged)
  }
}


extern __shared__ char s_array[];


__global__ void FIMCuda(float* d_triMem,float* d_triMemOut, int* d_vertMem,
  int* d_vertMemOutside, float* d_edgeMem0,float* d_edgeMem1,float* d_edgeMem2,
  float* d_speed, int* d_BlockSizes, int* d_con, int* ActiveList, int nActiveBlock,
  int maxNumTotalFaces, int maxNumVert,/*int nIter,*/ float m_StopDistance)
{
  uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;


  float* s_triMem  = (float*)s_array;


  // retrieve actual block index from the active list
  uint block_idx = ActiveList[list_idx];
  int block_size = d_BlockSizes[block_idx];


  ////////////////////////////////////////initialize shared memory//////////////////////////////////////////


  uint tri_base_addr = block_idx*maxNumTotalFaces*TRIMEMLENGTH;
  uint vert_base_addr = block_idx*maxNumVert*VERTMEMLENGTH;
  uint edge_base_addr = block_idx*maxNumTotalFaces;

  uint tx = threadIdx.x;



  short*   s_vertMem = (short*)&s_triMem[maxNumTotalFaces*TRIMEMLENGTH];



#pragma unroll
  for(int i = 0; i< TRIMEMLENGTH; i++)
  {
    s_triMem[tx*TRIMEMLENGTH + i] = d_triMem[tri_base_addr + tx * TRIMEMLENGTH + i];
  }

  if(tx < maxNumVert)
#pragma unroll
    for(int i = 0; i< VERTMEMLENGTH; i++)
    {
      s_vertMem[tx*VERTMEMLENGTH + i] = (short)(d_vertMem[vert_base_addr+tx*VERTMEMLENGTH + i] - tri_base_addr);

    }

  __syncthreads();


  /////////////////////////////////////////done shared memory copy//////////////////////////////////////////////////


  float a,b, delta, cosA, lamda1, lamda2, TC1, TC2;
  float TAB, TA, TB, TC;
  int C;
  float LenAB, LenBC, LenAC, LenCD, LenAD;
  float EdgeTA, EdgeTB;
  //float tmpOldTC;

  float oldT, newT;
  float oldValues0;
  float oldValues1;
  float oldValues2;
  //float oldValues[3];

  float speedI;



  //speedI = 1.0f;
  speedI = d_speed[edge_base_addr + tx];
  float edge0 = d_edgeMem0[edge_base_addr + tx];
  float edge1 = d_edgeMem1[edge_base_addr + tx];
  float edge2 = d_edgeMem2[edge_base_addr + tx];

#pragma unroll
  for(int iter=0; iter</*nIter*/NITER; iter++)
  {
    //
    // compute new value and unrolled the three computation
    //


    if(tx < block_size)
    {
      oldT = s_triMem[s_vertMem[tx*VERTMEMLENGTH + 0]];
    }

    __syncthreads();

    oldValues0 = s_triMem[tx*TRIMEMLENGTH + 0];
    oldValues1 = s_triMem[tx*TRIMEMLENGTH + 1];
    oldValues2 = s_triMem[tx*TRIMEMLENGTH + 2];

    if( !(oldValues0 >= LARGENUM && oldValues1 >= LARGENUM && oldValues2 >= LARGENUM) )  //if all values are large,break
    {



      C =  0;

      TA = oldValues1;
      TB = oldValues2;
      TC = oldValues0;




      TC1 = LARGENUM;
      TC2 = LARGENUM;




      TAB = TB - TA;

      //LenAB = s_triMem[tx*TRIMEMLENGTH + 1];     // the index should be tx*TRIMEMLENGTH + (i+1)%3
      //LenBC = s_triMem[tx*TRIMEMLENGTH + 2];     // the index should be tx*TRIMEMLENGTH + (i+2)%3
      //LenAC = s_triMem[tx*TRIMEMLENGTH + 0];     // the index should be tx*TRIMEMLENGTH + i

      LenAB = edge1;
      LenBC = edge2;
      LenAC = edge0;


      a = (speedI*speedI*LenAB*LenAB - TAB * TAB)*LenAB*LenAB;

      EdgeTA = TA + LenAC * speedI;
      EdgeTB = TB + LenBC * speedI;

      if (a > 0.0f)
      {

        cosA = (LenAC * LenAC + LenAB * LenAB - LenBC * LenBC) / (2.0f * LenAC * LenAB);

        b = 2.0f * LenAB * LenAC * cosA * (TAB * TAB - speedI*speedI*LenAB*LenAB);
        delta = 4.0f * LenAC * LenAC  * a *  TAB * TAB * (1.0f - cosA * cosA);

        lamda1 = (-b + sqrtf(delta))/(2.0f*a);
        lamda2 = (-b - sqrtf(delta))/(2.0f*a);

        if (lamda1>=0.0f &&lamda1<=1.0f)
        {
          LenAD = lamda1*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2*LenAC*LenAD*cosA);

          TC1 = lamda1*TAB+TA+LenCD*speedI;

        }
        if(lamda2>=0.0f && lamda2<=1.0f)
        {
          LenAD = lamda2*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2*LenAC*LenAD*cosA);

          TC2 = lamda2*TAB+TA+LenCD*speedI;

        }


        TC = MIN(TC, MIN(TC2, MIN(TC1,MIN(EdgeTA,EdgeTB))) );

      }

      else
      {
        TC = MIN(TC, MIN(EdgeTA,EdgeTB) );

      }

      s_triMem[tx*TRIMEMLENGTH + C] = TC;


      C = 1;

      TA = oldValues2;
      TB = oldValues0;
      TC = oldValues1;

      TC1 = LARGENUM;
      TC2 = LARGENUM;




      TAB = TB - TA;

      //LenAB = s_triMem[tx*TRIMEMLENGTH + A];     // the index should be tx*TRIMEMLENGTH + (i+1)%3
      //LenBC = s_triMem[tx*TRIMEMLENGTH + B];     // the index should be tx*TRIMEMLENGTH + (i+2)%3
      //LenAC = s_triMem[tx*TRIMEMLENGTH + C];     // the index should be tx*TRIMEMLENGTH + i

      LenAB = edge2;
      LenBC = edge0;
      LenAC = edge1;


      a = (speedI*speedI*LenAB*LenAB - TAB * TAB)*LenAB*LenAB;

      EdgeTA = TA + LenAC * speedI;
      EdgeTB = TB + LenBC * speedI;



      if (a > 0.0f)
      {

        cosA = (LenAC * LenAC + LenAB * LenAB - LenBC * LenBC) / (2.0f * LenAC * LenAB);
        b = 2.0f * LenAB * LenAC * cosA * (TAB * TAB - speedI*speedI*LenAB*LenAB);
        delta = 4.0f * LenAC * LenAC  * a *  TAB * TAB * (1.0f - cosA * cosA);

        lamda1 = (-b + sqrtf(delta))/(2.0f*a);
        lamda2 = (-b - sqrtf(delta))/(2.0f*a);

        if (lamda1>=0.0f &&lamda1<=1.0f)
        {
          LenAD = lamda1*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2*LenAC*LenAD*cosA);

          TC1 = lamda1*TAB+TA+LenCD*speedI;

        }
        if(lamda2>=0.0f &&lamda2<=1.0f)
        {
          LenAD = lamda2*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2.0f *LenAC*LenAD*cosA);

          TC2 = lamda2*TAB+TA+LenCD*speedI;

        }


        TC = MIN(TC, MIN(TC2, MIN(TC1,MIN(EdgeTA,EdgeTB))) );

      }

      else
      {
        TC = MIN(TC, MIN(EdgeTA,EdgeTB) );

      }

      s_triMem[tx*TRIMEMLENGTH + C] = TC;




      C = 2;
      TA = oldValues0;
      TB = oldValues1;
      TC = oldValues2;

      TC1 = LARGENUM;
      TC2 = LARGENUM;




      TAB = TB - TA;

      //LenAB = s_triMem[tx*TRIMEMLENGTH + A];     // the index should be tx*TRIMEMLENGTH + (i+1)%3
      //LenBC = s_triMem[tx*TRIMEMLENGTH + B];     // the index should be tx*TRIMEMLENGTH + (i+2)%3
      //LenAC = s_triMem[tx*TRIMEMLENGTH + C];     // the index should be tx*TRIMEMLENGTH + i

      LenAB = edge0;
      LenBC = edge1;
      LenAC = edge2;


      a = (speedI*speedI*LenAB*LenAB - TAB * TAB)*LenAB*LenAB;

      EdgeTA = TA + LenAC * speedI;
      EdgeTB = TB + LenBC * speedI;



      if (a > 0)
      {

        cosA = (LenAC * LenAC + LenAB * LenAB - LenBC * LenBC) / (2.0f * LenAC * LenAB);

        b = 2.0f * LenAB * LenAC * cosA * (TAB * TAB - speedI*speedI*LenAB*LenAB);
        delta = 4.0f * LenAC * LenAC  * a *  TAB * TAB * (1.0f - cosA * cosA);

        lamda1 = (-b + sqrtf(delta))/(2.0f*a);
        lamda2 = (-b - sqrtf(delta))/(2.0f*a);

        if (lamda1>=0.0f &&lamda1<=1.0f)
        {
          LenAD = lamda1*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2*LenAC*LenAD*cosA);

          TC1 = lamda1*TAB+TA+LenCD*speedI;

        }
        if(lamda2>=0.0f&&lamda2<=1.0f)
        {
          LenAD = lamda2*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2.0f *LenAC*LenAD*cosA);

          TC2 = lamda2*TAB+TA+LenCD*speedI;

        }


        TC = MIN(TC, MIN(TC2, MIN(TC1,MIN(EdgeTA,EdgeTB))) );

      }

      else
      {
        TC = MIN(TC, MIN(EdgeTA,EdgeTB) );

      }

      s_triMem[tx*TRIMEMLENGTH + C] = TC;

    }


    __syncthreads();


    TC = LARGENUM;


    if(tx < block_size)    //block_size is the vertices in this block and it is about warp size so there is no severe divergence
    {

      int tmp2;
      for(int j = 0; (j < VERTMEMLENGTH) && ( (tmp2 = s_vertMem[tx * VERTMEMLENGTH + j]) > -1); j++)                        // find the min and keep the old T for convergence check
      {

        TC = MIN(TC,s_triMem[tmp2]);

      }

      newT = TC;

      for(int j =0; (j < VERTMEMLENGTH) && ( (tmp2 = s_vertMem[tx * VERTMEMLENGTH + j]) > -1 ); j++) // update all the old to the min
      {
        s_triMem[tmp2] = TC;
      }
    }

    __syncthreads();
    ///////////////////////////////////////////////////////////////////////////////////////////////////

  }


  if(tx < block_size)
  {
    float residue = oldT - newT;
    if (residue < 0.) residue *= -1.;
    int tmpindex;
    d_con[block_idx*REDUCTIONSHARESIZE + tx] = (residue <= EPS) ? 1 : 0;

    // update gloal memory inside all the old to the min
    for(int j = 0; (j < VERTMEMLENGTH) &&
      ((tmpindex = s_vertMem[tx * VERTMEMLENGTH + j]) > -1); j++)  {
      d_triMemOut[tmpindex + tri_base_addr] = newT;
    }
    // update gloal memory outside all the old to the min
    for(int j = 0; (j < VERTMEMLENGTHOUTSIDE) &&
      ((tmpindex = d_vertMemOutside[block_idx*VERTMEMLENGTHOUTSIDE*maxNumVert + 
      tx * VERTMEMLENGTHOUTSIDE + j]) > -1 ); j++) {
      d_triMemOut[tmpindex] = newT;
    }
  }
}

extern __shared__ float s_run_check_neghbor_array[];



__global__ void run_check_neighbor(float* d_triMem,float* d_triMemOut, int* d_vertMem,
  int* d_vertMemOutside,float* d_edgeMem0,float* d_edgeMem1,float* d_edgeMem2, 
  float* d_speed,int* d_BlockSizes, int* d_con,int* d_ActiveList, int numOldActive ,
  int maxNumTotalFaces, int maxNumVert,int nTotalActive, int m_StopDistance)
{

  uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;


  if(list_idx < nTotalActive)
  {
    // retrieve actual block index from the active list
    uint block_idx = d_ActiveList[list_idx];
    int block_size = d_BlockSizes[block_idx];
    uint tri_base_addr = block_idx*maxNumTotalFaces*TRIMEMLENGTH;
    uint edge_base_addr = block_idx*maxNumTotalFaces;


    uint tx = threadIdx.x;


    uint vert_base_addr = block_idx*maxNumVert*VERTMEMLENGTH;




    float* s_triMem  = (float*)s_run_check_neghbor_array;
    short*   s_vertMem = (short*)&s_triMem[maxNumTotalFaces*TRIMEMLENGTH];

    float speedI;
    //speedI = 1.0f;
    speedI = d_speed[edge_base_addr + tx];




    float edge0 = d_edgeMem0[edge_base_addr + tx];
    float edge1 = d_edgeMem1[edge_base_addr + tx];
    float edge2 = d_edgeMem2[edge_base_addr + tx];

#pragma unroll
    for(int i = 0; i< TRIMEMLENGTH; i++)
    {
      s_triMem[tx*TRIMEMLENGTH + i] = d_triMem[tri_base_addr + tx * TRIMEMLENGTH + i];
    }

    if(tx < maxNumVert)
#pragma unroll
      for(int i = 0; i< VERTMEMLENGTH; i++)
      {
        s_vertMem[tx*VERTMEMLENGTH + i] = (short)(d_vertMem[vert_base_addr+tx*VERTMEMLENGTH + i] - tri_base_addr);

      }

    __syncthreads();


    /////////////////////////////////////////done shared memory copy//////////////////////////////////////////////////


    float a,b, delta, cosA, lamda1, lamda2, TC1, TC2;
    float TAB, TA, TB, TC;
    int C;
    float LenAB, LenBC, LenAC, LenCD, LenAD, EdgeTA, EdgeTB;

    float oldT = LARGENUM;
    float newT = LARGENUM;
    float oldValues0;
    float oldValues1;
    float oldValues2;


    //
    // compute new value
    //





    if(tx < block_size)
    {
      oldT = s_triMem[s_vertMem[tx*VERTMEMLENGTH + 0]];
    }

    oldValues0 = s_triMem[tx*TRIMEMLENGTH + 0];
    oldValues1 = s_triMem[tx*TRIMEMLENGTH + 1];
    oldValues2 = s_triMem[tx*TRIMEMLENGTH + 2];


    if( !(oldValues0 >= LARGENUM && oldValues1 >= LARGENUM && oldValues2 >= LARGENUM) )  //if all values are large,break
    {

      C =  0;

      TA = oldValues1;
      TB = oldValues2;
      TC = oldValues0;





      TC1 = LARGENUM;
      TC2 = LARGENUM;




      TAB = TB - TA;

      //LenAB = s_triMem[tx*TRIMEMLENGTH + A];     // the index should be tx*TRIMEMLENGTH + (i+1)%3
      //LenBC = s_triMem[tx*TRIMEMLENGTH + B];     // the index should be tx*TRIMEMLENGTH + (i+2)%3
      //LenAC = s_triMem[tx*TRIMEMLENGTH + C];     // the index should be tx*TRIMEMLENGTH + i

      LenAB = edge1;
      LenBC = edge2;
      LenAC = edge0;



      a = (speedI*speedI*LenAB*LenAB - TAB * TAB)*LenAB*LenAB;

      EdgeTA = TA + LenAC * speedI;
      EdgeTB = TB + LenBC * speedI;



      if (a > 0.0f)
      {

        cosA = (LenAC * LenAC + LenAB * LenAB - LenBC * LenBC) / (2.0f * LenAC * LenAB);
        b = 2.0f * LenAB * LenAC * cosA * (TAB * TAB - speedI*speedI*LenAB*LenAB);
        delta = 4.0f * LenAC * LenAC  * a *  TAB * TAB * (1.0f - cosA * cosA);

        lamda1 = (-b + sqrtf(delta))/(2.0f*a);
        lamda2 = (-b - sqrtf(delta))/(2.0f*a);

        if (lamda1>=0.0f&&lamda1<=1.0f)
        {
          LenAD = lamda1*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2.0f*LenAC*LenAD*cosA);

          TC1 = lamda1*TAB+TA+LenCD*speedI;

        }
        if(lamda2>=0.0f&&lamda2<=1.0f)
        {
          LenAD = lamda2*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2.0f*LenAC*LenAD*cosA);

          TC2 = lamda2*TAB+TA+LenCD*speedI;

        }


        TC = MIN(TC, MIN(TC2, MIN(TC1,MIN(EdgeTA,EdgeTB))) );

      }

      else
      {
        TC = MIN(TC, MIN(EdgeTA,EdgeTB) );

      }

      s_triMem[tx*TRIMEMLENGTH + C] = TC;


      C = 1;
      TA = oldValues2;
      TB = oldValues0;
      TC = oldValues1;


      TC1 = LARGENUM;
      TC2 = LARGENUM;




      TAB = TB - TA;

      //LenAB = s_triMem[tx*TRIMEMLENGTH + A];     // the index should be tx*TRIMEMLENGTH + (i+1)%3
      //LenBC = s_triMem[tx*TRIMEMLENGTH + B];     // the index should be tx*TRIMEMLENGTH + (i+2)%3
      //LenAC = s_triMem[tx*TRIMEMLENGTH + C];     // the index should be tx*TRIMEMLENGTH + i

      LenAB = edge2;
      LenBC = edge0;
      LenAC = edge1;


      a = (speedI*speedI*LenAB*LenAB - TAB * TAB)*LenAB*LenAB;

      EdgeTA = TA + LenAC * speedI;
      EdgeTB = TB + LenBC * speedI;



      if (a > 0.0f)
      {

        cosA = (LenAC * LenAC + LenAB * LenAB - LenBC * LenBC) / (2.0f * LenAC * LenAB);
        b = 2.0f * LenAB * LenAC * cosA * (TAB * TAB - speedI*speedI*LenAB*LenAB);
        delta = 4.0f * LenAC * LenAC  * a *  TAB * TAB * (1.0f - cosA * cosA);

        lamda1 = (-b + sqrtf(delta))/(2.0f *a);
        lamda2 = (-b - sqrtf(delta))/(2.0f *a);

        if (lamda1>=0.0f&&lamda1<=1.0f)
        {
          LenAD = lamda1*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2*LenAC*LenAD*cosA);

          TC1 = lamda1*TAB+TA+LenCD*speedI;

        }
        if(lamda2>=0.0f&&lamda2<=1.0f)
        {
          LenAD = lamda2*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2.0f *LenAC*LenAD*cosA);

          TC2 = lamda2*TAB+TA+LenCD*speedI;

        }


        TC = MIN(TC, MIN(TC2, MIN(TC1,MIN(EdgeTA,EdgeTB))) );

      }

      else
      {
        TC = MIN(TC, MIN(EdgeTA,EdgeTB) );

      }

      s_triMem[tx*TRIMEMLENGTH + C] = TC;




      C = 2;
      TA = oldValues0;
      TB = oldValues1;
      TC = oldValues2;


      TC1 = LARGENUM;
      TC2 = LARGENUM;




      TAB = TB - TA;

      //LenAB = s_triMem[tx*TRIMEMLENGTH + A];     // the index should be tx*TRIMEMLENGTH + (i+1)%3
      //LenBC = s_triMem[tx*TRIMEMLENGTH + B];     // the index should be tx*TRIMEMLENGTH + (i+2)%3
      //LenAC = s_triMem[tx*TRIMEMLENGTH + C];     // the index should be tx*TRIMEMLENGTH + i

      LenAB = edge0;
      LenBC = edge1;
      LenAC = edge2;


      a = (speedI*speedI*LenAB*LenAB - TAB * TAB)*LenAB*LenAB;

      EdgeTA = TA + LenAC * speedI;
      EdgeTB = TB + LenBC * speedI;



      if (a > 0.0f)
      {

        cosA = (LenAC * LenAC + LenAB * LenAB - LenBC * LenBC) / (2.0f * LenAC * LenAB);

        b = 2.0f * LenAB * LenAC * cosA * (TAB * TAB - speedI*speedI*LenAB*LenAB);
        delta = 4 * LenAC * LenAC  * a *  TAB * TAB * (1 - cosA * cosA);

        lamda1 = (-b + sqrtf(delta))/(2.0f*a);
        lamda2 = (-b - sqrtf(delta))/(2.0f*a);

        if (lamda1>=0&&lamda1<=1)
        {
          LenAD = lamda1*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2*LenAC*LenAD*cosA);

          TC1 = lamda1*TAB+TA+LenCD*speedI;

        }
        if(lamda2>=0.0f&&lamda2<=1.0f)
        {
          LenAD = lamda2*LenAB;

          LenCD = sqrtf(LenAC*LenAC+LenAD*LenAD-2.0f*LenAC*LenAD*cosA);

          TC2 = lamda2*TAB+TA+LenCD*speedI;

        }


        TC = MIN(TC, MIN(TC2, MIN(TC1,MIN(EdgeTA,EdgeTB))) );

      }

      else
      {
        TC = MIN(TC, MIN(EdgeTA,EdgeTB) );

      }

      s_triMem[tx*TRIMEMLENGTH + C] = TC;
    }


    __syncthreads();

    ///////////////////////////////////////////////////////////////////////////////////////////////////

    TC = LARGENUM;
    //int tmpcon = 1;

    if(tx < block_size)    //block_size is the vertices in this block and it is about warp size so there is no severe divergence
    {
      int tmp2;
      for(int j = 0; (j < VERTMEMLENGTH) && ( (tmp2 = s_vertMem[tx * VERTMEMLENGTH + j]) > -1); j++)                        // find the min and keep the old T for convergence check
      {

        TC = MIN(TC,s_triMem[tmp2/*+3*/]);

      }

      newT = TC;

      for(int j =0; (j < VERTMEMLENGTH) && ( (tmp2 = s_vertMem[tx * VERTMEMLENGTH + j]) > -1 ); j++) // update all the old to the min
      {
        s_triMem[tmp2] = TC;
      }
    }







    __syncthreads();



    if(tx < block_size)
    {
      float residue = oldT - newT;
      int tmpindex;
      d_con[block_idx*REDUCTIONSHARESIZE + tx] = (residue < EPS) ? 1 : 0;

      for(int j = 0; (j < VERTMEMLENGTH) && ( (tmpindex = s_vertMem[tx * VERTMEMLENGTH + j]) > -1 ); j++) // update gloal memory inside all the old to the min
      {

        d_triMemOut[tmpindex + tri_base_addr] = newT;
      }

      for(int j = 0; (j < VERTMEMLENGTHOUTSIDE) && ( (tmpindex = d_vertMemOutside[block_idx*VERTMEMLENGTHOUTSIDE*maxNumVert +  tx * VERTMEMLENGTHOUTSIDE + j]) > -1 ); j++) // update gloal memory outside all the old to the min
      {

        d_triMemOut[tmpindex] = newT;
      }

    }
    else if(tx < REDUCTIONSHARESIZE)
    {
      d_con[block_idx*REDUCTIONSHARESIZE + tx] = 1;   //assign the rest 1 so that in reduction, these values can be & with the effective values

    }


  }

}
