
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          modified to use only 7 floats for triMem
//1. #define TRIMEMLENGTH   7
//2. in FIMCuda and run_neighbor_check, add initilize old at the begining of iteration
//3. in FIMCuda and run_neighbor_check, s_triMem[tx*TRIMEMLENGTH + 3 + C] = TC after each iteration instead of s_triMem[tx*TRIMEMLENGTH + 6 + C] = TC
//4. in FIMCuda and run_neighbor_check, in the reconcile step, there should be no +3 in fetching the location of triMem
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meshFIM2d.h"
#include "Vec.h"
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cutil.h>
#include "CUDADefines.h"

#include <sstream>
#include <time.h>
#ifdef WIN32
#include <io.h>
#define unlink _unlink
#else
#include <unistd.h>
#endif
extern "C" {
#include <metis.h>
}



/////declaration for cuda kernels///////////////////////////
extern __global__ void run_reduction(int *con, int *blockCon,int* ActiveList, int nActiveBlock, int* blockSizes);
extern __global__ void FIMCuda(float* d_triMem,float* d_triMemOut, int* d_vertMem, int* d_vertMemOutside, float* d_edgeMem0,float* d_edgeMem1,float* d_edgeMem2,float* d_speed, int* d_BlockSizes, int* d_con, int* ActiveList, int nActiveBlock,int maxNumTotalFaces, int maxNumVert,/*int nIter, */float m_StopDistance);
extern __global__ void run_check_neighbor(float* d_triMem,float* d_triMemOut, int* d_vertMem,int* d_vertMemOutside,float* d_edgeMem0,float* d_edgeMem1,float* d_edgeMem2, float* d_speed, int* d_BlockSizes, int* d_con,int* d_ActiveList, int numOldActive ,int maxNumTotalFaces, int maxNumVert,int nTotalActive, int m_StopDistance);







#if __DEVICE_EMULATION__

bool InitCUDA(void){return true;}

#else
bool InitCUDA(void)
{
  int count = 0;
  int i = 0;

  cudaGetDeviceCount(&count);
  if(count == 0) {
    fprintf(stderr, "There is no device.\n");
    return false;
  }

  for(i = 0; i < count; i++) {
    cudaDeviceProp prop;
    if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
      if(prop.major >= 1) {
        break;
      }
    }
  }
  if(i == count) {
    fprintf(stderr, "There is no device supporting CUDA.\n");
    return false;
  }
  cudaSetDevice(i);

  printf("CUDA initialized.\n");
  return true;
}

#endif

/////////////////////////////////////////////////////////////////////////////

void meshFIM2d::GraphPartition_METIS(char* partfilename, int numBlock)  //read a metis result .mesh.npart.N file and store into PartitionLabel
{
  int numVert = m_meshPtr->vertices.size();
  //m_PartitionLabel = new int[numVert];
  m_PartitionLabel.resize(numVert);
  FILE* partFile = fopen(partfilename, "r+");
  if(partFile == NULL)
  {
    printf("NO part file found: %s\n", partfilename);
    exit(1);
  }


  for(int i = 0; i < numVert; i++)
  {
    fscanf(partFile, "%d", &m_PartitionLabel[i]);


  }

  srand( (unsigned)time( NULL ) );




  m_BlockSizes.resize(numBlock);

  //m_BlockSizes = new int[numBlock];
  for(int i =0; i<numBlock; i++)
    m_BlockSizes[i] = 0;

  float r,g,b;



  vector< Color > colors;
  colors.resize(numBlock);
  for(int i = 0; i< numBlock; i++)
  {
    r = rand()/(double)RAND_MAX;
    g = rand()/(double)RAND_MAX;
    b = rand()/(double)RAND_MAX;
    colors[i] = Color(r,g,b);
  }
  m_meshPtr->colors.resize(numVert);
  m_PartitionVerts.resize(numBlock);

  for(int i = 0; i<numVert; i++)
  {
    m_PartitionVerts[m_PartitionLabel[i]].push_back(i);


    m_BlockSizes[m_PartitionLabel[i]]++;
    m_meshPtr->colors[i] = colors[m_PartitionLabel[i]];

  }

  m_maxNumVert = 0;

  for(int i = 0 ; i < numBlock; i++)
  {

    m_maxNumVert = MAX(m_maxNumVert, m_BlockSizes[i]);
  }





}
//create .mesh file from trimesh faces and call partnmesh method 
//to partition and create intermediate mesh.npart.N file and then read this file
void meshFIM2d::GraphPartition_METIS2(int& numBlock, int maxNumBlockVerts, bool verbose)
{

  FILE * outf;

  outf = fopen("tmp.mesh", "w+");
  if(outf == NULL)
  {
    printf("Cannot open mesh file to write!!!!\n");
    exit(1);
  }
  int sz = m_meshPtr->faces.size();
  fprintf(outf,"%d 1\n", sz);

  for (int i=0;i<sz;i++)
    fprintf(outf, "%d %d %d\n",m_meshPtr->faces[i].v[0]+1,m_meshPtr->faces[i].v[1]+1,m_meshPtr->faces[i].v[2]+1);
  fclose(outf);

  int numVert = m_meshPtr->vertices.size();

  m_PartitionLabel.resize(numVert);

  char outputFileName[512];

  char meshfile[] = "tmp.mesh";

  if(numBlock == 0)
  {
    numBlock = numVert / maxNumBlockVerts;


    do{
      numBlock++;

      m_BlockSizes.resize(numBlock);
      for(int i=0; i< numBlock;i++)
      {
        m_BlockSizes[i] = 0;
      }
      partnmesh(meshfile,numBlock,verbose?1:0);

      sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);


      FILE* partFile = fopen(outputFileName, "r+");
      if(partFile == NULL)
      {
        printf("NO part file found: %s\n",outputFileName);
        exit(1);
      }

      for(int i = 0; i < numVert; i++)
      {
        fscanf(partFile, "%d", &m_PartitionLabel[i]);


      }

      for(int i = 0; i<numVert; i++)
      {


        m_BlockSizes[m_PartitionLabel[i]]++;

      }
      m_maxNumVert = 0;

      for(int i = 0 ; i < numBlock; i++)
      {

        m_maxNumVert = MAX(m_maxNumVert, m_BlockSizes[i]);
      }

      fclose(partFile);


      sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);
      unlink(outputFileName);
      sprintf(outputFileName, "tmp.mesh.epart.%d", numBlock);
      unlink(outputFileName);

    }while(m_maxNumVert != maxNumBlockVerts);
  }
  else
  {
    m_BlockSizes.resize(numBlock);
    for(int i=0; i< numBlock;i++)
    {
      m_BlockSizes[i] = 0;
    }

    partnmesh(meshfile,numBlock,verbose?1:0);

    sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);

    FILE* partFile = fopen(outputFileName, "r+");
    if(partFile == NULL)
    {
      printf("NO part file found: %s\n",outputFileName);
      exit(1);
    }

    for(int i = 0; i < numVert; i++)
    {
      fscanf(partFile, "%d", &m_PartitionLabel[i]);


    }

    for(int i = 0; i<numVert; i++)
    {


      m_BlockSizes[m_PartitionLabel[i]]++;

    }
    m_maxNumVert = 0;

    for(int i = 0 ; i < numBlock; i++)
    {

      m_maxNumVert = MAX(m_maxNumVert, m_BlockSizes[i]);
    }

    printf("max num vert is : %d\n", m_maxNumVert);
    fclose(partFile);


    sprintf(outputFileName, "tmp.mesh.npart.%d", numBlock);
    unlink(outputFileName);
    sprintf(outputFileName, "tmp.mesh.epart.%d", numBlock);
    unlink(outputFileName);

  }
  srand( (unsigned)time( NULL ) );
  printf("numBlock is : %d\n", numBlock);
  float r,g,b;
  vector< Color > colors;
  colors.resize(numBlock);
  for(int i = 0; i< numBlock; i++)
  {
    r = rand()/(double)RAND_MAX;
    g = rand()/(double)RAND_MAX;
    b = rand()/(double)RAND_MAX;
    colors[i] = Color(r,g,b);
  }
  m_meshPtr->colors.resize(numVert);
  m_PartitionVerts.resize(numBlock);

  for(int i = 0; i<numVert; i++)
  {
    m_PartitionVerts[m_PartitionLabel[i]].push_back(i);
    m_meshPtr->colors[i] = colors[m_PartitionLabel[i]];

  }
  unlink("tmp.mesh");
}

void meshFIM2d::GraphPartition_Square(int squareLength,int squareWidth, int blockLength, int blockWidth)
{
  int numVert = m_meshPtr->vertices.size();
  m_PartitionLabel.resize(numVert);

  int numBlockLength = (squareLength / blockLength);
  int numBlockWidth  = (squareWidth / blockWidth);
  int numBlock = numBlockLength * numBlockWidth;

  for(int i = 0; i< squareWidth; i++)
    for(int j =0; j< squareLength; j++)
    {
      m_PartitionLabel[i*squareLength+j] = (i/blockWidth) * numBlockLength + (j/blockLength);
    }

  m_BlockSizes.resize(numBlock);

  for(int i =0; i<numBlock; i++)
    m_BlockSizes[i] = 0;

  float r,g,b;
  vector< Color > colors;
  colors.resize(numBlock);
  for(int i = 0; i< numBlock; i++)
  {
    r = rand()/(double)RAND_MAX;
    g = rand()/(double)RAND_MAX;
    b = rand()/(double)RAND_MAX;
    colors[i] = Color(r,g,b);
  }
  m_meshPtr->colors.resize(numVert);
  m_PartitionVerts.resize(numBlock);

  for(int i = 0; i<numVert; i++)
  {
    m_PartitionVerts[m_PartitionLabel[i]].push_back(i);


    m_BlockSizes[m_PartitionLabel[i]]++;
    m_meshPtr->colors[i] = colors[m_PartitionLabel[i]];

  }

  m_maxNumVert = 0;

  for(int i = 0 ; i < numBlock; i++)
  {

    m_maxNumVert = MAX(m_maxNumVert, m_BlockSizes[i]);
  }
  printf("final number of blocks: %d\n", numBlock);
}

void meshFIM2d::PartitionFaces(int numBlock)
{
  ///////////////////step 3: partition faces//////////////////////////////////////
  m_PartitionFaces.resize(numBlock);
  m_PartitionNbFaces.resize(numBlock);

  int numFaces = m_meshPtr->faces.size();
  TriMesh::Face f;
  int labelv0;
  int labelv1;
  int labelv2;
  vector<TriMesh::Face> virtualfaces;
  vector<int> virtualFaceCnt;

  virtualFaceCnt.resize(numBlock);
  m_PartitionVirtualFaces.resize(numBlock);

  for(int i = 0; i< numBlock; i++)
    virtualFaceCnt[i] = 0;

  m_BlockNeighbor.resize(numBlock);

  for(int i = 0; i < numFaces; i++)
  {
    f = m_meshPtr->faces[i];
    int vfCnt = m_meshPtr->faceVirtualFaces[i].size();

    for(int k = 0 ; k < 3; k++)
    {
      if(!m_meshPtr->IsNonObtuse(f[k], f))
      {
        virtualFaceCnt[m_PartitionLabel[f[k]]] += vfCnt;
        m_PartitionVirtualFaces[m_PartitionLabel[f[k]]].insert(m_PartitionVirtualFaces[m_PartitionLabel[f[k]]].end(), m_meshPtr->faceVirtualFaces[i].begin(), m_meshPtr->faceVirtualFaces[i].end());
      }

    }

    labelv0 = m_PartitionLabel[f[0]];
    labelv1 = m_PartitionLabel[f[1]];
    labelv2 = m_PartitionLabel[f[2]];

    if(labelv0 == labelv1 && labelv1 == labelv2)
    {
      m_PartitionFaces[labelv0].push_back(i);
    }
    else if(labelv0 == labelv1 && labelv1 != labelv2)
    {
      m_PartitionNbFaces[labelv0].push_back(i);
      m_PartitionNbFaces[labelv2].push_back(i);

      m_BlockNeighbor[labelv0].insert(m_BlockNeighbor[labelv0].end(), labelv2);
      m_BlockNeighbor[labelv2].insert(m_BlockNeighbor[labelv2].end(), labelv0);
    }
    else if(labelv0 != labelv1 && labelv1 == labelv2)
    {
      m_PartitionNbFaces[labelv0].push_back(i);
      m_PartitionNbFaces[labelv2].push_back(i);

      m_BlockNeighbor[labelv0].insert(m_BlockNeighbor[labelv0].end(), labelv2);
      m_BlockNeighbor[labelv2].insert(m_BlockNeighbor[labelv2].end(), labelv0);
    }

    else if(labelv0 == labelv2 && labelv1 != labelv2)
    {
      m_PartitionNbFaces[labelv0].push_back(i);
      m_PartitionNbFaces[labelv1].push_back(i);

      m_BlockNeighbor[labelv0].insert(m_BlockNeighbor[labelv0].end(), labelv1);
      m_BlockNeighbor[labelv1].insert(m_BlockNeighbor[labelv1].end(), labelv0);
    }

    else      //all different
    {
      m_PartitionNbFaces[labelv0].push_back(i);
      m_PartitionNbFaces[labelv1].push_back(i);
      m_PartitionNbFaces[labelv2].push_back(i);

      m_BlockNeighbor[labelv0].insert(m_BlockNeighbor[labelv0].end(), labelv2);
      m_BlockNeighbor[labelv2].insert(m_BlockNeighbor[labelv2].end(), labelv0);
      m_BlockNeighbor[labelv0].insert(m_BlockNeighbor[labelv0].end(), labelv1);
      m_BlockNeighbor[labelv1].insert(m_BlockNeighbor[labelv1].end(), labelv0);
      m_BlockNeighbor[labelv1].insert(m_BlockNeighbor[labelv1].end(), labelv2);
      m_BlockNeighbor[labelv2].insert(m_BlockNeighbor[labelv2].end(), labelv1);
    }

  }

  vector<int> PartitionToltalFaces;
  PartitionToltalFaces.resize(numBlock);
  m_maxNumTotalFaces = 0;
  for(int j = 0; j < numBlock; j++)
  {


    PartitionToltalFaces[j] = m_PartitionFaces[j].size() + m_PartitionNbFaces[j].size() + virtualFaceCnt[j];
    m_maxNumTotalFaces = MAX(PartitionToltalFaces[j],m_maxNumTotalFaces );
  }
}

void meshFIM2d::GenerateData(int numBlock)
{
  int numVert = m_meshPtr->vertices.size();
  int numFaces=m_meshPtr->faces.size();

  if(!InitCUDA()) {
    exit(1);
  }

  index       *d_ActiveList= 0;
  int        *d_con;

  int*   d_con_forComputaion;
  int        *d_blockCon;
  float       *d_triMem;
  float      *d_edgeMem0;
  float      *d_edgeMem1;
  float      *d_edgeMem2;
  float*      d_speed;
  float      *d_triMemOut;
  int         *d_vertMem;
  int         *d_BlockSizes;

  index       *h_ActiveList= 0;    //list of active blocks
  int         *h_BlockLabel = 0;   //block active or not
  float       *h_triMem;
  float      *h_edgeMem0;
  float      *h_edgeMem1;
  float      *h_edgeMem2;
  float*      h_speed;
  int         *h_vertMem;
  int         *h_blockCon;

  int         *h_BlockSizes;

  /////////////////////////////malloc cpu memories///////////////////////////
  h_BlockLabel = (int*) malloc(sizeof(int) * numBlock);

  h_edgeMem0 = (float*)malloc(sizeof(float)  * m_maxNumTotalFaces * numBlock);
  h_edgeMem1 = (float*)malloc(sizeof(float)  * m_maxNumTotalFaces * numBlock);
  h_edgeMem2 = (float*)malloc(sizeof(float)  * m_maxNumTotalFaces * numBlock);
  h_speed    = (float*)malloc(sizeof(float)  * m_maxNumTotalFaces * numBlock);

  h_triMem = (float*)malloc(sizeof(float) * TRIMEMLENGTH * m_maxNumTotalFaces * numBlock);
  h_vertMem = (int*)malloc(sizeof(int) * VERTMEMLENGTH * m_maxNumVert * numBlock);
  h_BlockSizes = (int*)malloc(sizeof(int) * numBlock);
  h_blockCon = (int*)malloc(sizeof(int) * numBlock);
  /////////////////////////malloc gpu memories//////////////////////////////

  cudaSafeCall( cudaMalloc((void**) &d_con, sizeof(int) * numBlock * REDUCTIONSHARESIZE));

  cudaSafeCall( cudaMalloc((void**) &d_con_forComputaion, sizeof(int) * numBlock * REDUCTIONSHARESIZE));

  cudaSafeCall( cudaMalloc((void**) &d_blockCon,  sizeof(int) * numBlock));

  cudaSafeCall( cudaMalloc((void**) &d_triMem,  sizeof(float) * TRIMEMLENGTH * m_maxNumTotalFaces * numBlock));
  cudaSafeCall( cudaMalloc((void**) &d_triMemOut,  sizeof(float) * TRIMEMLENGTH * m_maxNumTotalFaces * numBlock));
  cudaSafeCall( cudaMalloc((void**) &d_edgeMem0,  sizeof(float)  * m_maxNumTotalFaces * numBlock));
  cudaSafeCall( cudaMalloc((void**) &d_edgeMem1,  sizeof(float)  * m_maxNumTotalFaces * numBlock));
  cudaSafeCall( cudaMalloc((void**) &d_edgeMem2,  sizeof(float)  * m_maxNumTotalFaces * numBlock));

  cudaSafeCall( cudaMalloc((void**) &d_speed,  sizeof(float)  * m_maxNumTotalFaces * numBlock));

  cudaSafeCall( cudaMalloc((void**) &d_vertMem, sizeof(int) * VERTMEMLENGTH * m_maxNumVert * numBlock));

  cudaSafeCall( cudaMalloc((void**) &d_BlockSizes, sizeof(int) * numBlock));
  /////////////////initialize cpu memories//////////////////////////////
  vector< vector<int> > blockVertMapping;
  blockVertMapping.resize(numVert);     //for each vertex, store the addresses where it appears in the global triMem array.
  for( int i = 0; i <  numBlock; i++)
  {
    int blockIdx = i * m_maxNumTotalFaces * TRIMEMLENGTH;
    int numPF = m_PartitionFaces[i].size();
    for(int j = 0; j< numPF; j++)
    {
      h_edgeMem0[i * m_maxNumTotalFaces + j] = m_meshPtr->faces[m_PartitionFaces[i][j]].edgeLens[0];
      h_edgeMem1[i * m_maxNumTotalFaces + j] = m_meshPtr->faces[m_PartitionFaces[i][j]].edgeLens[1];
      h_edgeMem2[i * m_maxNumTotalFaces + j] = m_meshPtr->faces[m_PartitionFaces[i][j]].edgeLens[2];

      h_triMem[blockIdx + j*TRIMEMLENGTH + 0] = LARGENUM;
      h_triMem[blockIdx + j*TRIMEMLENGTH + 1] = LARGENUM;
      h_triMem[blockIdx + j*TRIMEMLENGTH + 2] = LARGENUM;

      h_speed[i * m_maxNumTotalFaces + j]  =  m_meshPtr->faces[m_PartitionFaces[i][j]].speedInv;

      blockVertMapping[m_meshPtr->faces[m_PartitionFaces[i][j]][0]].push_back(blockIdx + j*TRIMEMLENGTH + 0);
      blockVertMapping[m_meshPtr->faces[m_PartitionFaces[i][j]][1]].push_back(blockIdx + j*TRIMEMLENGTH + 1);
      blockVertMapping[m_meshPtr->faces[m_PartitionFaces[i][j]][2]].push_back(blockIdx + j*TRIMEMLENGTH + 2);
    }
  }
  for( int i = 0; i <  numBlock; i++)
  {
    h_blockCon[i] = 1;

    h_BlockLabel[i] = m_BlockLabel[i];
    h_BlockSizes[i] = m_BlockSizes[i];
    int blockIdx = i * m_maxNumTotalFaces * TRIMEMLENGTH;

    int numPF = m_PartitionFaces[i].size();
    int numPNF = m_PartitionNbFaces[i].size();
    int numPVF = m_PartitionVirtualFaces[i].size();

    int k = 0;
    int l = 0;

    for(int j = numPF; j< m_maxNumTotalFaces; j++)
    {
      if( j < numPF + numPNF)
      {
        h_edgeMem0[i * m_maxNumTotalFaces + j] = m_meshPtr->faces[m_PartitionNbFaces[i][k]].edgeLens[0];
        h_edgeMem1[i * m_maxNumTotalFaces + j] = m_meshPtr->faces[m_PartitionNbFaces[i][k]].edgeLens[1];
        h_edgeMem2[i * m_maxNumTotalFaces + j]= m_meshPtr->faces[m_PartitionNbFaces[i][k]].edgeLens[2];

        h_triMem[blockIdx + j*TRIMEMLENGTH + 0] = LARGENUM;
        h_triMem[blockIdx + j*TRIMEMLENGTH + 1] = LARGENUM;
        h_triMem[blockIdx + j*TRIMEMLENGTH + 2] = LARGENUM;
        h_speed[i * m_maxNumTotalFaces + j] = m_meshPtr->faces[m_PartitionNbFaces[i][k]].speedInv;

        blockVertMapping[m_meshPtr->faces[m_PartitionNbFaces[i][k]][0]].push_back(blockIdx + j*TRIMEMLENGTH + 0);
        blockVertMapping[m_meshPtr->faces[m_PartitionNbFaces[i][k]][1]].push_back(blockIdx + j*TRIMEMLENGTH + 1);
        blockVertMapping[m_meshPtr->faces[m_PartitionNbFaces[i][k]][2]].push_back(blockIdx + j*TRIMEMLENGTH + 2);

        k++;
      }
      else if (j < numPF + numPNF + numPVF)
      {
        h_edgeMem0[i * m_maxNumTotalFaces + j]= m_PartitionVirtualFaces[i][l].edgeLens[0];
        h_edgeMem1[i * m_maxNumTotalFaces + j]= m_PartitionVirtualFaces[i][l].edgeLens[1];
        h_edgeMem2[i * m_maxNumTotalFaces + j] = m_PartitionVirtualFaces[i][l].edgeLens[2];

        h_triMem[blockIdx + j*TRIMEMLENGTH + 0] = LARGENUM;
        h_triMem[blockIdx + j*TRIMEMLENGTH + 1] = LARGENUM;
        h_triMem[blockIdx + j*TRIMEMLENGTH + 2] = LARGENUM;
        h_speed[i * m_maxNumTotalFaces + j]  =m_PartitionVirtualFaces[i][l].speedInv;

        blockVertMapping[m_PartitionVirtualFaces[i][l][0]].push_back(blockIdx + j*TRIMEMLENGTH + 0);
        blockVertMapping[m_PartitionVirtualFaces[i][l][1]].push_back(blockIdx + j*TRIMEMLENGTH + 1);
        blockVertMapping[m_PartitionVirtualFaces[i][l][2]].push_back(blockIdx + j*TRIMEMLENGTH + 2);

        l++;
      }
      else
      {
        h_triMem[blockIdx + j*TRIMEMLENGTH + 0] = LARGENUM;
        h_triMem[blockIdx + j*TRIMEMLENGTH + 1] = LARGENUM;
        h_triMem[blockIdx + j*TRIMEMLENGTH + 2] = LARGENUM;
      }
    }
  }

  m_maxNumVertMapping = 0;
  for(int i =0; i < numVert; i++)
  {
    int blockIndex = m_PartitionLabel[i];
    int tmp = blockVertMapping[i][0];
    int maxi = (blockIndex+1) * m_maxNumTotalFaces * TRIMEMLENGTH;
    int mini = blockIndex * m_maxNumTotalFaces * TRIMEMLENGTH;
    if(  ( tmp< mini) || (tmp >= maxi) )
    {
      for(int j =0; j < blockVertMapping[i].size(); j++)
        if(blockVertMapping[i][j] >= mini && blockVertMapping[i][j] < maxi )
        {
          int swaptmp = tmp;
          blockVertMapping[i][0] = blockVertMapping[i][j];
          blockVertMapping[i][j] = swaptmp;
          break;

        }
    }
    m_maxNumVertMapping = MAX(m_maxNumVertMapping, blockVertMapping[i].size());
  }

  for(int i =0; i < numVert; i++)
  {
    int blockIndex = m_PartitionLabel[i];
    int tmp = blockVertMapping[i][0];
    int maxi = (blockIndex+1) * m_maxNumTotalFaces * TRIMEMLENGTH;
    int mini = blockIndex * m_maxNumTotalFaces * TRIMEMLENGTH;
    if(  ( tmp< mini) || (tmp >= maxi) )
    {
      printf("beyond");
    }
  }

  vector< vector<int> > blockVertMappingInside;
  vector< vector<int> > blockVertMappingOutside;

  blockVertMappingInside.resize(numVert);
  blockVertMappingOutside.resize(numVert);

  for(int i = 0; i< numBlock; i++)
  {
    int triIdx =  i * TRIMEMLENGTH * m_maxNumTotalFaces;

    for(int m  = 0; m < m_PartitionVerts[i].size(); m++)
    {

      vector<int> tmp = blockVertMapping[m_PartitionVerts[i][m]];


      for(int n = 0; n < tmp.size(); n++)
      {
        if( tmp[n] >= triIdx + 0  && tmp[n] < triIdx + m_maxNumTotalFaces*TRIMEMLENGTH)
          blockVertMappingInside[m_PartitionVerts[i][m]].push_back(tmp[n]);
        else
        {
          blockVertMappingOutside[m_PartitionVerts[i][m]].push_back(tmp[n]);

        }
      }

    }
  }

  int maxVertMappingInside = 0;
  int maxVertMappingOutside = 0;
  for(int i =0; i< numVert; i++)
  {
    maxVertMappingInside = MAX(maxVertMappingInside, (blockVertMappingInside[i].size()));
    maxVertMappingOutside = MAX(maxVertMappingInside, (blockVertMappingOutside[i].size()));
  }

  printf("maxVertMappingInside is: %d\n",maxVertMappingInside);
  printf("maxVertMappingOutside is: %d\n",maxVertMappingOutside);


  for(int i = 0; i< numBlock; i++)
  {
    int vertIdx =  i * VERTMEMLENGTH * m_maxNumVert;

    for(int m  = 0; m < m_PartitionVerts[i].size(); m++)
    {

      int tmpsize = blockVertMappingInside[m_PartitionVerts[i][m]].size();

      int n = 0;
      for(; n < tmpsize; n++)
        h_vertMem[vertIdx + m*VERTMEMLENGTH + n] = blockVertMappingInside[m_PartitionVerts[i][m]][n];
      for(;n<VERTMEMLENGTH; n++)

        h_vertMem[vertIdx + m*VERTMEMLENGTH + n] = -1 + i*m_maxNumTotalFaces*TRIMEMLENGTH;

    }

    for(int m = m_PartitionVerts[i].size() * VERTMEMLENGTH; m < m_maxNumVert * VERTMEMLENGTH; m++)
    {
      h_vertMem[vertIdx + m] = -1 + i*m_maxNumTotalFaces*TRIMEMLENGTH;
    }
  }


  int* h_vertMemOutside = (int*)malloc(m_maxNumVert * numBlock * VERTMEMLENGTHOUTSIDE * sizeof(int));
  int* d_vertMemOutside;
  cudaSafeCall( cudaMalloc((void**) &d_vertMemOutside, m_maxNumVert * numBlock * VERTMEMLENGTHOUTSIDE * sizeof(int) ) );

  for(int i = 0; i< numBlock; i++)
  {
    int vertIdx =  i * VERTMEMLENGTHOUTSIDE * m_maxNumVert;

    for(int m  = 0; m < m_PartitionVerts[i].size(); m++)
    {

      int tmpsize = blockVertMappingOutside[m_PartitionVerts[i][m]].size();

      int n = 0;
      for(; n < tmpsize; n++)
        h_vertMemOutside[vertIdx + m*VERTMEMLENGTHOUTSIDE + n] = blockVertMappingOutside[m_PartitionVerts[i][m]][n];
      for(;n<VERTMEMLENGTHOUTSIDE; n++)
        h_vertMemOutside[vertIdx + m*VERTMEMLENGTHOUTSIDE + n] = -1;

    }

    for(int m = m_PartitionVerts[i].size() * VERTMEMLENGTHOUTSIDE; m < m_maxNumVert * VERTMEMLENGTHOUTSIDE; m++)
    {
      h_vertMemOutside[vertIdx + m] = -1;
    }
  }

  h_ActiveList = (int*)malloc(sizeof(int)*numBlock);
  cudaSafeCall( cudaMalloc((void**) &d_ActiveList, sizeof(int) * numBlock));

  //////////////////////////////////////////////////////////////////////////////////
  vector<int>  nb;
  int numActive;

  for( int i = 0; i <  numBlock; i++)
  {

    h_blockCon[i] = 1;

    h_BlockLabel[i] = m_BlockLabel[i];
    h_BlockSizes[i] = m_BlockSizes[i];
  }

  //////////////initialize the seed points for h_triMem////////////////////////////////////

  for(int i = 0; i< m_SeedPoints.size(); i++)
  {
    int seed = m_SeedPoints[i];
    int seedBelongToBlock = m_PartitionLabel[seed];
    h_blockCon[seedBelongToBlock] = 0;
    for(int j = 0; j < blockVertMapping[seed].size(); j++)
    {
      h_triMem[blockVertMapping[seed][j]] = 0.0;
    }
  }
  /////////////copy triMem and verMem to a vector just for debugging/////////////////
  vector<float> vec_triMem;
  vector<int>   vec_vertMem;
  vector<int>   vec_vertMemOutside;

  vec_triMem.resize(TRIMEMLENGTH * m_maxNumTotalFaces * numBlock);
  vec_vertMem.resize(VERTMEMLENGTH * m_maxNumVert * numBlock);
  vec_vertMemOutside.resize(VERTMEMLENGTHOUTSIDE * m_maxNumVert * numBlock);
  for(int i =0; i < TRIMEMLENGTH * m_maxNumTotalFaces * numBlock; i++)
    vec_triMem[i] = h_triMem[i];

  for(int i = 0; i< VERTMEMLENGTH * m_maxNumVert * numBlock; i++)
    vec_vertMem[i] = h_vertMem[i];
  for(int i = 0; i< VERTMEMLENGTHOUTSIDE * m_maxNumVert * numBlock; i++)
    vec_vertMemOutside[i] = h_vertMemOutside[i];
  ////////////////////////////////////////////////////////////////////////////

  cudaSafeCall( cudaMemcpy( d_triMem,h_triMem, sizeof(float) * m_maxNumTotalFaces * numBlock * TRIMEMLENGTH, cudaMemcpyHostToDevice));

  numActive =m_ActiveBlocks.size();

  set<int>::iterator activeiter = m_ActiveBlocks.begin();
  for(int i =0; activeiter !=  m_ActiveBlocks.end(); activeiter++)
    h_ActiveList[i++] = *activeiter;

  cudaEvent_t start, stop, startCopy, stopCopy;
  cudaEventCreate(&start);
  cudaEventCreate(&startCopy);
  cudaEventCreate(&stopCopy);
  cudaEventCreate(&stop);
  cudaEventRecord(startCopy,0);


  //////////////////copy to gpu memories///////////////////////////////

  cudaSafeCall( cudaMemcpy( d_triMem,h_triMem, sizeof(float) * m_maxNumTotalFaces * numBlock * TRIMEMLENGTH, cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_triMemOut,h_triMem, sizeof(float) * m_maxNumTotalFaces * numBlock * TRIMEMLENGTH, cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_edgeMem0,h_edgeMem0, sizeof(float) * m_maxNumTotalFaces * numBlock , cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_edgeMem1,h_edgeMem1, sizeof(float) * m_maxNumTotalFaces * numBlock , cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_edgeMem2,h_edgeMem2, sizeof(float) * m_maxNumTotalFaces * numBlock , cudaMemcpyHostToDevice));

  cudaSafeCall( cudaMemcpy( d_speed,h_speed, sizeof(float) * m_maxNumTotalFaces * numBlock , cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_vertMem,h_vertMem, sizeof(int) * m_maxNumVert * numBlock * VERTMEMLENGTH, cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_vertMemOutside,h_vertMemOutside, sizeof(int) * m_maxNumVert * numBlock * VERTMEMLENGTHOUTSIDE, cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_BlockSizes,h_BlockSizes, sizeof(int) * numBlock, cudaMemcpyHostToDevice));
  cudaSafeCall( cudaMemcpy( d_blockCon,h_blockCon, sizeof(int) * numBlock, cudaMemcpyHostToDevice));

  printf("max number of triangles per block: %d\n", m_maxNumTotalFaces);
  int nTotalIter = 0;

  cudaEventRecord(start,0);

  int totalIterationNumber = 0;

  while ( numActive > 0)
  {
    ///////////step 1: run solver /////////////////////////////////////////////////////////////

    nTotalIter++;

    totalIterationNumber += numActive;
    printf("number of active block: %d\n", numActive);

    dim3 dimGrid(numActive, 1);
    dim3 dimBlock(m_maxNumTotalFaces, 1);

    cudaSafeCall( cudaMemcpy( d_ActiveList,h_ActiveList,sizeof(int) * numBlock, cudaMemcpyHostToDevice));

    FIMCuda<<<dimGrid, dimBlock, m_maxNumTotalFaces*TRIMEMLENGTH*sizeof(float)+m_maxNumVert*VERTMEMLENGTH*sizeof(short)>>>( d_triMem,d_triMemOut, d_vertMem,d_vertMemOutside,d_edgeMem0,d_edgeMem1,d_edgeMem2, d_speed, d_BlockSizes, d_con,d_ActiveList, numActive,m_maxNumTotalFaces, m_maxNumVert, m_StopDistance);
    cudaCheckErrors();

    //////////////////////step 2: reduction////////////////////////////////////////////////

    dimBlock = dim3(REDUCTIONSHARESIZE / 2 , 1);
    run_reduction<<<dimGrid, dimBlock/*, sizeof(int)*m_maxNumVert*/>>>(d_con, d_blockCon,d_ActiveList, numActive, d_BlockSizes);
    cudaCheckErrors();

    //////////////////////////////////////////////////////////////////
    // 3. check neighbor tiles of converged tile
    // Add any active block of neighbor of converged block is inserted
    // to the list

    cudaSafeCall( cudaMemcpy(h_blockCon, d_blockCon, numBlock*sizeof(int), cudaMemcpyDeviceToHost) );

    int nOldActiveBlock = numActive;

    for(uint i=0; i<nOldActiveBlock; i++)
    {
      // check neighbors of current active tile
      uint currBlkIdx = h_ActiveList[i];

      if(h_blockCon[currBlkIdx]) // not active : converged
      {
        //h_BlockLabel[currBlkIdx] == FARP;
        set<int> nb = m_BlockNeighbor[currBlkIdx];

        set<int>::iterator iter;
        for( iter = nb.begin(); iter != nb.end() ; iter++)
        {
          int currIdx = *iter;

          if(h_BlockLabel[currIdx] == FARP)
          {
            h_BlockLabel[currIdx] = ACTIVE;
            h_ActiveList[numActive++] = currIdx;
            //m_ActiveBlocks.insert(m_ActiveBlocks.end(), currIdx);
          }
        }
      }

    }
    //////////////////////////////////////////////////////////////////
    // 4. run solver only once for neighbor blocks of converged block
    // current active list contains active blocks and neighbor blocks of
    // any converged blocks
    //

    cudaSafeCall( cudaMemcpy(d_ActiveList, h_ActiveList, numActive*sizeof(int), cudaMemcpyHostToDevice) );

    dimGrid = dim3(numActive, 1);
    dimBlock = dim3(m_maxNumTotalFaces, 1);


    run_check_neighbor<<< dimGrid, dimBlock, m_maxNumTotalFaces*TRIMEMLENGTH*sizeof(float)+m_maxNumVert*VERTMEMLENGTH*sizeof(short)>>>(d_triMemOut, d_triMem,d_vertMem,d_vertMemOutside, d_edgeMem0,d_edgeMem1,d_edgeMem2,d_speed , d_BlockSizes, d_con,d_ActiveList, nOldActiveBlock ,m_maxNumTotalFaces, m_maxNumVert,numActive, m_StopDistance);
    cudaCheckErrors();


    //////////////////////////////////////////////////////////////////
    // 5. reduction

    dimGrid = dim3(numActive, 1);
    dimBlock = dim3(REDUCTIONSHARESIZE / 2 , 1);

    run_reduction<<<dimGrid, dimBlock/*, sizeof(int)*m_maxNumVert*/>>>(d_con, d_blockCon,d_ActiveList,numActive, d_BlockSizes);
    cudaCheckErrors();

    //////////////////////////////////////////////////////////////////
    // 6. update active list
    // read back active volume from the device and add
    // active block to active list on the host memory

    numActive = 0;

    cudaSafeCall( cudaMemcpy(h_blockCon, d_blockCon, numBlock*sizeof(int), cudaMemcpyDeviceToHost) );
    for(uint i=0; i<numBlock; i++)
    {
      if(!h_blockCon[i]) // false : activate block (not converged)
      {
        h_BlockLabel[i] = ACTIVE;
        h_ActiveList[numActive++] = i;
      }
      else h_BlockLabel[i] = FARP;
    }
    ////////////////////////DEBUG see if this is possible
    cudaSafeCall( cudaMemcpy(h_triMem, d_triMem,sizeof(float) *
          m_maxNumTotalFaces * numBlock * TRIMEMLENGTH , cudaMemcpyDeviceToHost) );
    for(int i =0; i < numVert; i++) {
      m_meshPtr->vertT[i] =  h_triMem[blockVertMapping[i][0]];
    }
    std::stringstream fname;
    fname << "result2D-";
    size_t digits = 3;
    int tmp = nTotalIter;

    while (tmp / 10 > 0) {
      digits --;
      tmp /=10;
    }
    for (size_t i = 0; i < digits; i++)
      fname << "0";
    fname << nTotalIter << ".txt";
    FILE * resultfile = fopen(fname.str().c_str(), "w+");
    for(int i = 0; i < numVert; i++)
    {
      fprintf(resultfile, "%.8f\n", m_meshPtr->vertT[i]);
    }
    fclose(resultfile);
    ////////////////////////////////END DEBUG
  }

  cudaSafeCall( cudaThreadSynchronize() );

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);

  cudaSafeCall( cudaMemcpy(h_triMem, d_triMem,sizeof(float) * m_maxNumTotalFaces * numBlock * TRIMEMLENGTH , cudaMemcpyDeviceToHost) );

  cudaSafeCall( cudaThreadSynchronize() );

  cudaEventRecord(stopCopy,0);
  cudaEventSynchronize(stopCopy);

  float totalTime, totalAndCopyTime;
  cudaEventElapsedTime(&totalTime, start, stop);
  cudaEventElapsedTime(&totalAndCopyTime, startCopy, stopCopy);


  cudaCheckErrors();


  printf("Total Processing time: %f (ms)\n", totalTime);
  printf("Total Processing time and copy time: %f (ms)\n", totalAndCopyTime);
  printf("The iteration number: %d\n", nTotalIter );
  printf("The total iteration number: %d\n", totalIterationNumber );
  printf("The total localsolver calls per vertex: %f\n", totalIterationNumber*m_maxNumTotalFaces*(NITER+1)*3.0 / (float)numVert);

  vec_triMem.resize(m_maxNumTotalFaces * numBlock * 3);
  float maxVertT = 0;
  for(int i = 0 ; i <  m_maxNumTotalFaces * numBlock; i++)
  {

    vec_triMem[3*i + 0] =  h_triMem[i*TRIMEMLENGTH + 0];
    vec_triMem[3*i + 1] =  h_triMem[i*TRIMEMLENGTH + 1];
    vec_triMem[3*i + 2] =  h_triMem[i*TRIMEMLENGTH + 2];

    if(h_triMem[i*TRIMEMLENGTH + 0] >= LARGENUM)
      vec_triMem[3*i + 0] = -2;
    if(h_triMem[i*TRIMEMLENGTH + 1] >= LARGENUM)
      vec_triMem[3*i + 1] = -2;
    if(h_triMem[i*TRIMEMLENGTH + 2] >= LARGENUM)
      vec_triMem[3*i + 2] = -2;


    maxVertT = MAX(maxVertT,MAX(vec_triMem[3*i + 2] , MAX(vec_triMem[3*i + 1] , vec_triMem[3*i + 0])));
  }

  int vertIndex = 0;

  for(int i =0; i < numVert; i++)
  {
    m_meshPtr->vertT[i] =  h_triMem[blockVertMapping[i][0]];
    if(m_meshPtr->vertT[i] == maxVertT)
      vertIndex = i;
  }
  printf("The maximun vertT is: %f, the vert index is: %d \n", maxVertT,vertIndex );
  cudaSafeCall( cudaFree(d_ActiveList));
  cudaSafeCall( cudaFree(d_triMem));
  cudaSafeCall( cudaFree(d_vertMem));
  cudaSafeCall( cudaFree(d_edgeMem0));
  cudaSafeCall( cudaFree(d_edgeMem1));
  cudaSafeCall( cudaFree(d_edgeMem2));
  cudaSafeCall( cudaFree(d_speed));
  cudaSafeCall( cudaFree(d_con));
  cudaSafeCall( cudaFree(d_blockCon));
  free(h_ActiveList);
  free(h_edgeMem0);
  free(h_edgeMem1);
  free(h_edgeMem2);
  free(h_speed);
  free(h_triMem);
  free(h_vertMem);
  free(h_BlockLabel);
  free(h_blockCon);
  free(h_BlockSizes);
}
