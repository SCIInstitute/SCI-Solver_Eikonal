
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          modified to use only 7 floats for triMem
//1. #define TRIMEMLENGTH   7
//2. in FIMCuda and run_neighbor_check, add initilize old at the begining of iteration
//3. in FIMCuda and run_neighbor_check, s_triMem[tx*TRIMEMLENGTH + 3 + C] = TC after each iteration instead of s_triMem[tx*TRIMEMLENGTH + 6 + C] = TC
//4. in FIMCuda and run_neighbor_check, in the reconcile step, there should be no +3 in fetching the location of triMem
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meshFIM3d.h"
#include "Vec.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef WIN32
#include <io.h>
#define unlink _unlink
#else
#include <unistd.h>
#endif
#include "CUDADefines.h"
#include <time.h>
#include <mycutil.h>
extern "C" {
#include <metis.h>
}


/////declaration for cuda kernels///////////////////////////
extern __global__ void run_reduction(bool *con, bool *blockCon, int* ActiveList, int nActiveBlock, int* blockSizes);


extern __global__ void FIMCuda(float3* d_tetMem0, float3* d_tetMem1, float4* d_tetT, float* d_vertT, float* d_speedInv, int* d_vertMem, int* d_vertMemOutside,
    int* d_BlockSizes, bool* d_con, int* d_ActiveList,
    int m_maxNumInVert, int m_maxVertMappingInside, int m_maxNumOutVertMapping, int nIter);

extern __global__ void CopyOutBack(float4* d_tetT, float* d_vertT, int* d_vertMem, int* d_vertMemOutside, int* d_BlockSizes, int* d_ActiveList, int m_maxNumInVert, int m_maxNumTotalTets, int m_maxVertMappingInside, int m_maxVertMappingOutside);


extern __global__ void run_check_neghbor(float3* d_tetMem0, float3* d_tetMem1, float4* d_tetT, float* d_speedInv, int* d_vertMem, int* d_vertMemOutside,
    int* d_BlockSizes, bool* d_con, int* d_ActiveList,
    int m_maxNumInVert, int m_maxVertMappingInside, int m_maxNumOutVertMapping);

#if __DEVICE_EMULATION__

bool InitCUDA(void)
{
  return true;
}


#else

bool InitCUDA(void)
{
  int count = 0;
  int i = 0;

  cudaGetDeviceCount(&count);
  if(count == 0)
  {
    fprintf(stderr, "There is no device.\n");
    return false;
  }

  for(i = 0; i < count; i++)
  {
    cudaDeviceProp prop;
    if(cudaGetDeviceProperties(&prop, i) == cudaSuccess)
    {
      if(prop.major >= 1)
      {
        break;
      }
    }
  }
  if(i == count)
  {
    fprintf(stderr, "There is no device supporting CUDA.\n");
    return false;
  }

  cudaDeviceProp props;
  cudaSafeCall(cudaSetDevice(0));

  cudaSafeCall(cudaGetDeviceProperties(&props, 0));

  printf("Device 0: \"%s\" with Compute %d.%d capability\n",  props.name, props.major, props.minor);

  printf("CUDA initialized.\n");
  return true;
}

#endif

/////////////////////////////////////////////////////////////////////////////
//create .mesh file from trimesh faces and call partnmesh function
//to partition and create intermediate mesh.npart.N file and then read this file
void meshFIM3d::GraphPartition_METIS2(int& numBlock, int maxNumBlockVerts, bool verbose)
{

  FILE * outf;

  outf = fopen("tmp.mesh", "w+");
  if(outf == NULL)
  {
    printf("Cannot open mesh file to write!!!!\n");
    exit(1);
  }
  int sz = m_meshPtr->tets.size();

  fprintf(outf, "%d 2\n", sz);

  for(int i = 0; i < sz; i++)
    fprintf(outf, "%d %d %d %d\n", m_meshPtr->tets[i].v[0] + 1, m_meshPtr->tets[i].v[1] + 1, m_meshPtr->tets[i].v[2] + 1, m_meshPtr->tets[i].v[3] + 1);

  fclose(outf);


  int numVert = m_meshPtr->vertices.size();

  m_PartitionLabel.resize(numVert);

  char outputFileName[512];

  char meshfile[] = "tmp.mesh";

  if(numBlock == 0)
  {
    numBlock = MAX(numVert / maxNumBlockVerts - 10, 1);


    do
    {
      numBlock++;

      m_BlockSizes.resize(numBlock);
      for(int i = 0; i < numBlock; i++)
      {
        m_BlockSizes[i] = 0;
      }
      partnmesh(meshfile,numBlock,verbose?1:0);

      sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);


      FILE* partFile = fopen(outputFileName, "r+");
      if(partFile == NULL)
      {
        printf("NO part file found!!!!\n");
        exit(1);
      }
      for(int i = 0; i < numVert; i++)
      {
        fscanf(partFile, "%d", &m_PartitionLabel[i]);
      }

      for(int i = 0; i < numVert; i++)
      {
        m_BlockSizes[m_PartitionLabel[i]]++;
      }
      m_maxNumInVert = 0;

      for(int i = 0; i < numBlock; i++)
      {

        m_maxNumInVert = MAX(m_maxNumInVert, m_BlockSizes[i]);
      }

      fclose(partFile);

      sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);
      unlink(outputFileName);
      sprintf(outputFileName, "tmp.mesh.epart.%d", numBlock);
      unlink(outputFileName);

    }
    while(m_maxNumInVert != maxNumBlockVerts);
  }
  else
  {
    m_BlockSizes.resize(numBlock);
    for(int i = 0; i < numBlock; i++)
    {
      m_BlockSizes[i] = 0;
    }

    partnmesh(meshfile, numBlock,verbose?1:0);

    sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);

    FILE* partFile = fopen(outputFileName, "r+");
    if(partFile == NULL)
    {
      printf("NO part file found!!!!\n");
      exit(1);
    }

    for(int i = 0; i < numVert; i++)
    {
      fscanf(partFile, "%d", &m_PartitionLabel[i]);
    }

    for(int i = 0; i < numVert; i++)
    {
      m_BlockSizes[m_PartitionLabel[i]]++;
    }
    m_maxNumInVert = 0;

    for(int i = 0; i < numBlock; i++)
    {
      m_maxNumInVert = MAX(m_maxNumInVert, m_BlockSizes[i]);
    }

    printf("max num vert is : %d\n", m_maxNumInVert);
    fclose(partFile);

    sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);
    unlink(outputFileName);
    sprintf(outputFileName, "tmp.mesh.epart.%d", numBlock);
    unlink(outputFileName);
  }

  srand((unsigned)time(NULL));

  printf("numBlock is : %d\n", numBlock);

  m_PartitionInVerts.resize(numBlock);

  for(int i = 0; i < numVert; i++)
  {
    m_PartitionInVerts[m_PartitionLabel[i]].push_back(i);
  }
  unlink("tmp.mesh");
}

void meshFIM3d::GraphPartition_Square(int squareLength, int squareWidth, int squareHeight, int blockLength, int blockWidth, int blockHeight)
{
  int numVert = m_meshPtr->vertices.size();
  m_PartitionLabel.resize(numVert);

  int numBlockLength = (squareLength / blockLength);
  int numBlockWidth = (squareWidth / blockWidth);
  int numBlockHeight = (squareHeight / blockHeight);
  int numBlock = numBlockLength * numBlockWidth*numBlockHeight;


  for(int k = 0; k < squareHeight; k++)
    for(int i = 0; i < squareWidth; i++)
      for(int j = 0; j < squareLength; j++)
      {

        int index = k * squareLength * squareWidth + i * squareLength + j;
        int k2 = k;
        int i2 = i;
        int j2 = j;
        m_PartitionLabel[index] = (k2 / blockHeight) * numBlockLength * numBlockWidth + (i2 / blockWidth) * numBlockLength + (j2 / blockLength);
      }

  m_BlockSizes.resize(numBlock);

  for(int i = 0; i < numBlock; i++)
    m_BlockSizes[i] = 0;

  m_PartitionInVerts.resize(numBlock);

  for(int i = 0; i < numVert; i++)
  {
    m_PartitionInVerts[m_PartitionLabel[i]].push_back(i);
    m_BlockSizes[m_PartitionLabel[i]]++;
  }

  m_maxNumInVert = 0;

  for(int i = 0; i < numBlock; i++)
  {
    m_maxNumInVert = MAX(m_maxNumInVert, m_BlockSizes[i]);
  }
  printf("final number of blocks: %d\n", numBlock);
}

void meshFIM3d::GenerateData(void)
{

  printf("Start GenerateData!\n");

  int numVert = m_meshPtr->vertices.size();

  if(!InitCUDA())
  {
    exit(1);
  }

  float* h_tetMem0;
  float* h_tetMem1;
  float* h_tetT;
  float* h_vertT;
  int* h_vertMem;
  int* h_vertMemOutside;
  bool* h_blockCon;
  int* h_BlockSizes;
  int* h_BlockLabel;
  vector<int> h_ActiveList;
  vector<int> h_ActiveListNew;

  int* d_ActiveList = 0;
  bool* d_con;
  bool* d_blockCon;
  float3* d_tetMem0;
  float3* d_tetMem1;
  float4* d_tetT;
  float* d_vertT;
  float* d_speedInv;
  int* d_vertMem;
  int* d_vertMemOutside;
  int* d_BlockSizes;

  GetTetMem(h_tetMem0, h_tetMem1, h_tetT);
  GetVertMem(h_vertMem, h_vertMemOutside);
  h_vertT = (float*)malloc(sizeof(float)* m_maxNumInVert * m_numBlock);


  h_blockCon = (bool*)malloc(sizeof(bool) * m_numBlock);
  h_BlockLabel = (int*)malloc(sizeof(int)* m_numBlock);
  h_BlockSizes = (int*)malloc(sizeof(int)* m_numBlock);
  memset(h_blockCon, 1, sizeof(bool) * m_numBlock);
  for(int i = 0; i < m_numBlock; i++) {
    h_BlockLabel[i] = FARP;
    h_BlockSizes[i] = m_BlockSizes[i];
  }

  ////////////////////initialize the seed points for h_tetT//////////////////////////

  printf("Seed size is %lu, source block is %d\n", m_SeedPoints.size(),
      m_PartitionLabel.empty()?-1:
      (m_PartitionLabel[m_SeedPoints.empty()?0:m_SeedPoints[0]]));
  for(int i = 0; i < m_SeedPoints.size(); i++)
  {
    int seed = m_SeedPoints[i];
    int seedBelongToBlock = m_PartitionLabel[seed];
    m_ActiveBlocks.insert(m_ActiveBlocks.end(), seedBelongToBlock);
    h_blockCon[seedBelongToBlock] = false;
    h_BlockLabel[seedBelongToBlock] = ACTIVE;
    // int blockIdx = seedBelongToBlock * m_maxNumTotalFaces * TRIMEMLENGTH;
    for(int j = 0; j < m_blockVertMapping[seed].size(); j++)
    {
      h_tetT[m_blockVertMapping[seed][j]] = 0.0;
    }
  }

  int numActive = m_ActiveBlocks.size();

  printf("Active block number is %d.\n", numActive);


  h_ActiveList.resize(m_numBlock);

  set<int>::iterator activeiter = m_ActiveBlocks.begin();
  for(int i = 0; activeiter != m_ActiveBlocks.end(); activeiter++)
    h_ActiveList[i++] = *activeiter;


  unsigned int timerstart, timerend = 0;


  ///////////////////////malloc GPU memory/////////////////////////////////
  cudaSafeCall((cudaMalloc((void**)&d_con, sizeof(bool) * m_numBlock * m_maxNumInVert)));
  cudaSafeCall((cudaMalloc((void**)&d_tetMem0, sizeof(float)* 3 * m_maxNumTotalTets * m_numBlock)));
  cudaSafeCall((cudaMalloc((void**)&d_tetMem1, sizeof(float)* 3 * m_maxNumTotalTets * m_numBlock)));
  cudaSafeCall((cudaMalloc((void**)&d_tetT, sizeof(float)* 4 * m_maxNumTotalTets * m_numBlock)));
  cudaSafeCall((cudaMalloc((void**)&d_vertT, sizeof(float)* m_maxNumInVert * m_numBlock)));
  cudaSafeCall( cudaMalloc( (void**) &d_speedInv, sizeof(float) * m_maxNumTotalTets * m_numBlock) );
  cudaSafeCall((cudaMalloc((void**)&d_vertMem, sizeof(int)* m_maxNumInVert * m_numBlock * m_maxVertMappingInside)));
  cudaSafeCall((cudaMalloc((void**)&d_vertMemOutside, sizeof(int)* m_maxNumInVert * m_numBlock * m_maxVertMappingOutside)));
  cudaSafeCall((cudaMalloc((void**)&d_BlockSizes, sizeof(int)* m_numBlock)));
  cudaSafeCall((cudaMalloc((void**)&d_blockCon, sizeof(bool) * m_numBlock)));
  cudaSafeCall((cudaMalloc((void**)&d_ActiveList, sizeof(int)* m_numBlock)));


  //////////////////copy to gpu memories///////////////////////////////
  cudaSafeCall((cudaMemcpy(d_tetMem0, h_tetMem0, sizeof(float)* 3 * m_maxNumTotalTets * m_numBlock, cudaMemcpyHostToDevice)));
  cudaSafeCall((cudaMemcpy(d_tetMem1, h_tetMem1, sizeof(float)* 3 * m_maxNumTotalTets * m_numBlock, cudaMemcpyHostToDevice)));
  cudaSafeCall((cudaMemcpy(d_tetT, h_tetT, sizeof(float)* 4 * m_maxNumTotalTets * m_numBlock, cudaMemcpyHostToDevice)));
  //  cudaSafeCall( cudaMemcpy( d_speedInv,h_speedInv, sizeof(float) * m_maxNumTotalTets * m_numBlock , cudaMemcpyHostToDevice));
  cudaSafeCall((cudaMemcpy(d_vertMem, h_vertMem, sizeof(int)* m_maxNumInVert * m_numBlock * m_maxVertMappingInside, cudaMemcpyHostToDevice)));
  cudaSafeCall((cudaMemcpy(d_vertMemOutside, h_vertMemOutside, sizeof(int)* m_maxNumInVert * m_numBlock * m_maxVertMappingOutside, cudaMemcpyHostToDevice)));
  cudaSafeCall((cudaMemcpy(d_BlockSizes, h_BlockSizes, sizeof(int)* m_numBlock, cudaMemcpyHostToDevice)));
  cudaSafeCall((cudaMemcpy(d_blockCon, h_blockCon, sizeof(bool) * m_numBlock, cudaMemcpyHostToDevice)));
  cudaSafeCall((cudaMemset(d_vertT, 0, sizeof(float)* m_maxNumInVert * m_numBlock)));

  int nTotalIter = 0;
  int nIter = 2;


  cudaFuncSetCacheConfig(FIMCuda, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(run_check_neghbor, cudaFuncCachePreferShared);

  vector< vector<float> > tmp_h_verrT;
  vector< vector<float> > tmp_h_verrT2;
  tmp_h_verrT.resize(m_numBlock);
  tmp_h_verrT2.resize(m_numBlock);

  int totalIterationNumber = 0;
  timerstart = clock();

  printf("Start Iteration!!\n");

  int maxActive = 0;
  while(numActive > 0)
  {
    maxActive = MAX(maxActive, numActive);
    printf("nTotalIter = %d, numActive = %d\n", nTotalIter, numActive);
    ///////step 1: run solver /////////////////////////////////////
    nTotalIter++;
    totalIterationNumber += numActive;
    dim3 dimGrid(numActive, 1);
    //dim3 dimBlock(520, 1);
    dim3 dimBlock(m_maxNumTotalTets, 1);
    cudaSafeCall(cudaMemcpy(d_ActiveList, &h_ActiveList[0], sizeof(int)* m_numBlock, cudaMemcpyHostToDevice));
    int sharedSize = sizeof(float)* 4 * m_maxNumTotalTets + sizeof(int)* m_maxNumInVert * m_maxVertMappingInside;
    cudaSafeCall((FIMCuda << <dimGrid, dimBlock, sharedSize >> >(d_tetMem0, d_tetMem1, d_tetT, d_vertT, d_speedInv, d_vertMem, d_vertMemOutside,
            d_BlockSizes, d_con, d_ActiveList, m_maxNumInVert, m_maxVertMappingInside, m_maxVertMappingOutside, nIter)));

    dimBlock = dim3(m_maxNumInVert, 1);
    CopyOutBack << <dimGrid, dimBlock >> >(d_tetT, d_vertT, d_vertMem, d_vertMemOutside, d_BlockSizes, d_ActiveList, m_maxNumInVert, m_maxNumTotalTets, m_maxVertMappingInside, m_maxVertMappingOutside);
    //////////////////////step 2: reduction////////////////////////////////////////////////
    dimBlock = dim3(m_maxNumInVert, 1);
    run_reduction << <dimGrid, dimBlock >> >(d_con, d_blockCon, d_ActiveList, numActive, d_BlockSizes);

    //////////////////////////////////////////////////////////////////
    // 3. check neighbor tiles of converged tile
    // Add any active block of neighbor of converged block is inserted
    // to the list
    cudaSafeCall(cudaMemcpy(h_blockCon, d_blockCon, m_numBlock * sizeof(bool), cudaMemcpyDeviceToHost));

    int nOldActiveBlock = numActive;

    for(uint i = 0; i < nOldActiveBlock; i++)
    {
      // check neighbors of current active tile
      uint currBlkIdx = h_ActiveList[i];

      if(h_blockCon[currBlkIdx]) // not active : converged
      {
        set<int> nb = m_BlockNeighbor[currBlkIdx];

        set<int>::iterator iter;
        for( iter = nb.begin(); iter != nb.end() ; iter++)
        {
          int currIdx = *iter;

          if(h_BlockLabel[currIdx] == FARP)
          {
            h_BlockLabel[currIdx] = ACTIVE;
            h_ActiveList[numActive++] = currIdx;
          }
        }
      }
    }
    //      //////////////////////////////////////////////////////////////////
    //      // 4. run solver only once for neighbor blocks of converged block
    //      // current active list contains active blocks and neighbor blocks of
    //      // any converged blocks

    cudaSafeCall(cudaMemcpy(d_ActiveList, &h_ActiveList[0], numActive * sizeof(int), cudaMemcpyHostToDevice));
    dimGrid = dim3(numActive, 1);
    dimBlock = dim3(m_maxNumTotalTets, 1);

    sharedSize = sizeof(float4) * m_maxNumTotalTets + sizeof(int)* m_maxNumInVert * m_maxVertMappingInside;

    run_check_neghbor << <dimGrid, dimBlock, sharedSize >> >(d_tetMem0, d_tetMem1, d_tetT, d_speedInv, d_vertMem, d_vertMemOutside,
        d_BlockSizes, d_con, d_ActiveList, m_maxNumInVert, m_maxVertMappingInside, m_maxVertMappingOutside);

    ////////////////////////////////////////////////////////////////
    // 5. reduction
    ///////////////////////////////////////////////////////////////
    dimGrid = dim3(numActive, 1);
    dimBlock = dim3(m_maxNumInVert, 1);
    run_reduction << <dimGrid, dimBlock >> >(d_con, d_blockCon, d_ActiveList, numActive, d_BlockSizes);

    //////////////////////////////////////////////////////////////////
    // 6. update active list
    // read back active volume from the device and add
    // active block to active list on the host memory

    numActive = 0;

    cudaSafeCall(cudaMemcpy(h_blockCon, d_blockCon, m_numBlock * sizeof(bool), cudaMemcpyDeviceToHost));
    for(uint i = 0; i < m_numBlock; i++)
    {
      if(!h_blockCon[i]) // false : activate block (not converged)
      {
        h_BlockLabel[i] = ACTIVE;
        h_ActiveList[numActive++] = i;
      }
      else h_BlockLabel[i] = FARP;
    }
  }

  cudaSafeCall(cudaThreadSynchronize());
  timerend = clock();
  double duration = (double)(timerend - timerstart) / CLOCKS_PER_SEC;

  printf("Computing time : %.10lf s\n",duration);

  cudaSafeCall(cudaMemcpy(h_vertT, d_vertT, sizeof(float)* m_maxNumInVert * m_numBlock, cudaMemcpyDeviceToHost));

  cudaSafeCall(cudaThreadSynchronize());

  printf("num of max active %d\n", maxActive);

  m_meshPtr->vertT.resize(1);
  m_meshPtr->vertT[0].resize(numVert);

#ifdef _DEBUG

  for(int i = 0; i < m_numBlock; i++)
  {

    tmp_h_verrT[i].resize(m_maxNumInVert);
    for(int j = 0; j < m_maxNumInVert; j++)
    {

      tmp_h_verrT[i][j] = h_vertT[i * m_maxNumInVert + j];


    }
  }

#endif

  //FILE* resultfile = fopen("resultXX.txt", "w+");
  for(int i = 0; i < m_numBlock; i++)
  {
    for(int j = 0; j < m_PartitionInVerts[i].size(); j++)
    {
      //fprintf(resultfile, "%.8f\n", h_vertT[i * m_maxNumInVert + j]);
      m_meshPtr->vertT[0][m_PartitionInVerts[i][j]] = h_vertT[i * m_maxNumInVert + j];
    }
  }
  //fclose(resultfile);

  FILE * resultfile = fopen("result.txt", "w+");
  for(int i = 0; i < numVert; i++)
  {
    fprintf(resultfile, "%.8f\n", m_meshPtr->vertT[0][i]);
  }

  fclose(resultfile);



  printf("The iteration number: %d\n", nTotalIter);
  printf("The total iteration number: %d\n", totalIterationNumber);

  cudaSafeCall(cudaFree(d_con));
  cudaSafeCall(cudaFree(d_blockCon));
  cudaSafeCall(cudaFree(d_BlockSizes));
  cudaSafeCall(cudaFree(d_speedInv));

  free(h_blockCon);
  free(h_BlockSizes);
}

void meshFIM3d::PartitionTets(int numBlock)
{
  ///////////////////////////////////step 3: partition faces//////////////////////////////////////
  printf("Start PartitionTets ...");
  m_PartitionTets.resize(numBlock);
  m_PartitionNbTets.resize(numBlock);

  int numTets = m_meshPtr->tets.size();
  int numVerts = m_meshPtr->vertices.size();
  TetMesh::Tet t;

  vector<TetMesh::Tet> virtualTets;
  vector<int> virtualTetCnt;

  virtualTetCnt.resize(numBlock);
  m_PartitionVirtualTets.resize(numBlock);
  set<int> labels;

  for(int i = 0; i < numTets; i++)
  {
    t = m_meshPtr->tets[i];
    int vfCnt = m_meshPtr->tetVirtualTets[i].size();


    int obtusevert = t.obtuseV;
    if(obtusevert >= 0)
    {
      virtualTetCnt[m_PartitionLabel[t[obtusevert]]] += vfCnt;
      m_PartitionVirtualTets[m_PartitionLabel[t[obtusevert]]].insert(m_PartitionVirtualTets[m_PartitionLabel[t[obtusevert]]].end(), m_meshPtr->tetVirtualTets[i].begin(), m_meshPtr->tetVirtualTets[i].end());
    }
    labels.clear();
    for(int m = 0; m < 4; m++)
      labels.insert(labels.begin(), m_PartitionLabel[t[m]]);
    if(labels.size() == 1)
    {
      m_PartitionTets[*(labels.begin())].push_back(i);
    }
    else if(labels.size() > 1)
    {
      set<int>::iterator it = labels.begin();
      for(set<int>::iterator it = labels.begin(); it != labels.end(); it++)
      {
        m_PartitionNbTets[*it].push_back(i);
      }
    }
    else
      printf("Error!!\n");
  }

  vector<int> PartitionToltalTets;
  PartitionToltalTets.resize(numBlock);
  m_maxNumTotalTets = 0;
  for(int j = 0; j < numBlock; j++)
  {
    PartitionToltalTets[j] = m_PartitionTets[j].size() + m_PartitionNbTets[j].size() + virtualTetCnt[j];
    m_maxNumTotalTets = MAX(PartitionToltalTets[j], m_maxNumTotalTets);
  }

  printf("m_maxNumTotalTets is %d\n", m_maxNumTotalTets);


  //calculate block neighbors.
  m_BlockNeighbor.resize(numBlock);
  for(int i = 0; i < numVerts; i++)
  {
    vector<int> nbs = m_meshPtr->neighbors[i];
    for(int j = 0; j < nbs.size(); j++)
    {
      int nb = nbs[j];
      if(m_PartitionLabel[nb] != m_PartitionLabel[i])
        m_BlockNeighbor[m_PartitionLabel[i]].insert(m_BlockNeighbor[m_PartitionLabel[i]].end(), m_PartitionLabel[nb]);
    }

  }
  printf("done!\n");
}

bool meshFIM3d::gettetmem(vector<float>& tetmem, TetMesh::Tet t)
{
  bool needswap = false;
  tetmem.resize(6);
  point A = m_meshPtr->vertices[t[0]];
  point B = m_meshPtr->vertices[t[1]];
  point C = m_meshPtr->vertices[t[2]];
  point D = m_meshPtr->vertices[t[3]];

  point AB = B - A;
  point AC = C - A;
  point AD = D - A;

  AC = C - A;
  AD = D - A;
  point BC = C - B;
  point CD = D - C;
  point BD = D - B;

  tetmem[0] = vMv(AC, t.M, BC);
  tetmem[1] = vMv(BC, t.M, CD);
  tetmem[2] = vMv(AC, t.M, CD);
  tetmem[3] = vMv(AD, t.M, BD);
  tetmem[4] = vMv(AC, t.M, AD);
  tetmem[5] = vMv(BC, t.M, BD);

  return needswap;

}

void meshFIM3d::GetTetMem(float* &h_tetMem0, float* &h_tetMem1, float* &h_tetT)
{
  h_tetMem0 = (float*)malloc(3 * sizeof(float)* m_maxNumTotalTets * m_numBlock);
  h_tetMem1 = (float*)malloc(3 * sizeof(float)* m_maxNumTotalTets * m_numBlock);
  h_tetT = (float*)malloc(4 * sizeof(float)* m_maxNumTotalTets * m_numBlock);

  int numTets = m_meshPtr->tets.size();

  int numVert = m_meshPtr->vertices.size();

  m_blockVertMapping.resize(numVert); //for each vertex, store the addresses where it appears in the global triMem array.

  TetMesh::Tet t;


  for(int i = 0; i < m_numBlock; i++)
  {
    int blockIdx = i * m_maxNumTotalTets * 3;
    int numPF = m_PartitionTets[i].size();
    for(int j = 0; j < numPF; j++)
    {

      t = m_meshPtr->tets[m_PartitionTets[i][j]];
      vector<float> tetmem;
      bool needswap = gettetmem(tetmem, t);

      h_tetMem0[blockIdx + j * 3 + 0] = tetmem[0];
      h_tetMem0[blockIdx + j * 3 + 1] = tetmem[1];
      h_tetMem0[blockIdx + j * 3 + 2] = tetmem[2];

      h_tetMem1[blockIdx + j * 3 + 0] = tetmem[3];
      h_tetMem1[blockIdx + j * 3 + 1] = tetmem[4];
      h_tetMem1[blockIdx + j * 3 + 2] = tetmem[5];

      h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 0] = LARGENUM_TET;
      h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 1] = LARGENUM_TET;
      h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 2] = LARGENUM_TET;
      h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 3] = LARGENUM_TET;

      m_blockVertMapping[t[0]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 0);
      m_blockVertMapping[t[3]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 3);

      if(needswap)
      {
        m_blockVertMapping[t[1]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 2);
        m_blockVertMapping[t[2]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 1);
      }
      else
      {
        m_blockVertMapping[t[1]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 1);
        m_blockVertMapping[t[2]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 2);

      }
    }

  }

  for(int i = 0; i < m_numBlock; i++)
  {
    int blockIdx = i * m_maxNumTotalTets * 3;

    int numPF = m_PartitionTets[i].size();
    int numPNF = m_PartitionNbTets[i].size();
    int numPVF = m_PartitionVirtualTets[i].size();

    int k = 0;
    int l = 0;

    for(int j = numPF; j < m_maxNumTotalTets; j++)
    {

      if(j < numPF + numPNF)
      {

        vector<float> tetmem;
        t = m_meshPtr->tets[m_PartitionNbTets[i][k]];
        bool needswap = gettetmem(tetmem, t);

        h_tetMem0[blockIdx + j * 3 + 0] = tetmem[0];
        h_tetMem0[blockIdx + j * 3 + 1] = tetmem[1];
        h_tetMem0[blockIdx + j * 3 + 2] = tetmem[2];

        h_tetMem1[blockIdx + j * 3 + 0] = tetmem[3];
        h_tetMem1[blockIdx + j * 3 + 1] = tetmem[4];
        h_tetMem1[blockIdx + j * 3 + 2] = tetmem[5];

        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 0] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 1] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 2] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 3] = LARGENUM_TET;

        m_blockVertMapping[t[0]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 0);
        m_blockVertMapping[t[3]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 3);
        if(needswap)
        {
          m_blockVertMapping[t[1]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 2);
          m_blockVertMapping[t[2]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 1);
        }
        else
        {
          m_blockVertMapping[t[1]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 1);
          m_blockVertMapping[t[2]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 2);
        }
        k++;
      }
      else if(j < numPF + numPNF + numPVF)
      {
        t = m_PartitionVirtualTets[i][l];
        vector<float> tetmem;
        bool needswap = gettetmem(tetmem, t);

        h_tetMem0[blockIdx + j * 3 + 0] = tetmem[0];
        h_tetMem0[blockIdx + j * 3 + 1] = tetmem[1];
        h_tetMem0[blockIdx + j * 3 + 2] = tetmem[2];

        h_tetMem1[blockIdx + j * 3 + 0] = tetmem[3];
        h_tetMem1[blockIdx + j * 3 + 1] = tetmem[4];
        h_tetMem1[blockIdx + j * 3 + 2] = tetmem[5];

        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 0] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 1] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 2] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 3] = LARGENUM_TET;

        m_blockVertMapping[t[0]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 0);
        m_blockVertMapping[t[3]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 3);
        if(needswap)
        {
          m_blockVertMapping[t[1]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 2);
          m_blockVertMapping[t[2]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 1);
        }
        else
        {
          m_blockVertMapping[t[1]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 1);
          m_blockVertMapping[t[2]].push_back(i * m_maxNumTotalTets * 4 + j * 4 + 2);
        }
        l++;
      }
      else
      {
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 0] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 1] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 2] = LARGENUM_TET;
        h_tetT[i * m_maxNumTotalTets * 4 + j * 4 + 3] = LARGENUM_TET;
      }
    }
  }
}

void meshFIM3d::GetVertMem(int* &h_vertMem, int* &h_vertMemOutside)
{

  int numVert = m_meshPtr->vertices.size();

  m_blockVertMappingInside.resize(numVert);
  m_blockVertMappingOutside.resize(numVert);

  m_maxNumVertMapping = 0;

  for(int i = 0; i < m_numBlock; i++)
  {
    int triIdx = i * TETMEMLENGTH * m_maxNumTotalTets;

    for(int m = 0; m < m_PartitionInVerts[i].size(); m++)
    {

      m_maxNumVertMapping = MAX(m_maxNumVertMapping, m_blockVertMapping[i].size());

      vector<int> tmp = m_blockVertMapping[m_PartitionInVerts[i][m]];


      for(int n = 0; n < tmp.size(); n++)
      {
        if(tmp[n] >= triIdx + 0 && tmp[n] < triIdx + m_maxNumTotalTets * TETMEMLENGTH)
          m_blockVertMappingInside[m_PartitionInVerts[i][m]].push_back(tmp[n]);
        else
        {
          m_blockVertMappingOutside[m_PartitionInVerts[i][m]].push_back(tmp[n]);
        }
      }
    }
  }

  m_maxVertMappingInside = 0;
  m_maxVertMappingOutside = 0;
  for(int i = 0; i < numVert; i++)
  {
    m_maxVertMappingInside = MAX(m_maxVertMappingInside, (m_blockVertMappingInside[i].size()));
    m_maxVertMappingOutside = MAX(m_maxVertMappingOutside, (m_blockVertMappingOutside[i].size()));
  }

  printf("maxVertMappingInside is: %d\n", m_maxVertMappingInside);
  printf("maxVertMappingOutside is: %d\n", m_maxVertMappingOutside);

  h_vertMem = (int*)malloc(sizeof(int)* m_maxVertMappingInside * m_maxNumInVert * m_numBlock);
  for(int i = 0; i < m_numBlock; i++)
  {
    int vertIdx = i * m_maxVertMappingInside * m_maxNumInVert;

    for(int m = 0; m < m_PartitionInVerts[i].size(); m++)
    {

      int tmpsize = m_blockVertMappingInside[m_PartitionInVerts[i][m]].size();

      int n = 0;
      for(; n < tmpsize; n++)
        h_vertMem[vertIdx + m * m_maxVertMappingInside + n] = m_blockVertMappingInside[m_PartitionInVerts[i][m]][n];
      for(; n < m_maxVertMappingInside; n++)
        h_vertMem[vertIdx + m * m_maxVertMappingInside + n] = -1 + i * m_maxNumTotalTets * TETMEMLENGTH;

    }

    for(int m = m_PartitionInVerts[i].size() * m_maxVertMappingInside; m < m_maxNumInVert * m_maxVertMappingInside; m++)
    {
      h_vertMem[vertIdx + m] = -1 + i * m_maxNumTotalTets*TETMEMLENGTH;
    }
  }


  h_vertMemOutside = (int*)malloc(m_maxNumInVert * m_numBlock * m_maxVertMappingOutside * sizeof(int));

  for(int i = 0; i < m_numBlock; i++)
  {
    int vertIdx = i * m_maxVertMappingOutside * m_maxNumInVert;

    for(int m = 0; m < m_PartitionInVerts[i].size(); m++)
    {

      int tmpsize = m_blockVertMappingOutside[m_PartitionInVerts[i][m]].size();

      int n = 0;
      for(; n < tmpsize; n++)
        h_vertMemOutside[vertIdx + m * m_maxVertMappingOutside + n] = m_blockVertMappingOutside[m_PartitionInVerts[i][m]][n];
      for(; n < m_maxVertMappingOutside; n++)
        h_vertMemOutside[vertIdx + m * m_maxVertMappingOutside + n] = -1;

    }

    for(int m = m_PartitionInVerts[i].size() * m_maxVertMappingOutside; m < m_maxNumInVert * m_maxVertMappingOutside; m++)
    {
      h_vertMemOutside[vertIdx + m] = -1;
    }
  }
  printf("Done GetVertMem!!\n");
}
