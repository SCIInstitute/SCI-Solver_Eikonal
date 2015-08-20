//-------------------------------------------------------------------
//
//  Copyright (C) 2015
//  Scientific Computing & Imaging Institute
//  University of Utah
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "meshFIM3d.h"
#include <math.h>
#include "tetgen.h"
#include "tetmesh.h"

#include <time.h>

using std::string;

int main(int argc, char *argv[])
{
  //Verbose option
  bool verbose = false;
  int numBlock = 10003;
  size_t maxIter = 100;
  //input filename (minus extension)
  std::string filename;
  for (int i = 0; i < argc; i++)
    if (strcmp(argv[i],"-v") == 0) {
      verbose = true;
    } else if (strcmp(argv[i],"-m") == 0) {
      if (i+1 >= argc) break;
      maxIter = atoi(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-b") == 0) {
      if (i+1 >= argc) break;
      numBlock = atoi(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-i") == 0) {
      if (i+1 >= argc) break;
      filename = std::string(argv[i+1]);
      i++;
    }
  if (filename.empty())
    filename = "example_data/CubeMesh_size16";
  clock_t starttime = clock(), endtime;

  tetgenio in, addin, bgmin,out;
  in.load_tetmesh((char*)filename.c_str());

  TetMesh themesh;
  themesh.init(
      in.pointlist,
      in.numberofpoints,
      in.trifacelist,
      in.numberoffacets,
      in.tetrahedronlist,
      in.numberoftetrahedra,
      in.numberoftetrahedronattributes,
      in.tetrahedronattributelist );

  themesh.need_neighbors();
  themesh.need_adjacenttets();
  themesh.need_tet_virtual_tets();

  meshFIM3d* FIMPtr = new meshFIM3d;

  vector<int> seedPointList(1,0/*133152*//*20181*//*2184*/);//20181 for unstruc_s5
  //seedPointList.push_back(10);
  //seedPointList.push_back(20);
  //seedPointList.push_back(30);
  //seedPointList.push_back(40);
  FIMPtr->SetSeedPoint(seedPointList);

  FIMPtr->SetMesh(&themesh);
  //FIMPtr->FindSeedPoint();
  FIMPtr->InitSpeedMat();

  int squareLength = 16;
  int squareWidth = 16;
  int squareDepth = 16;
  int squareBlockLength = 4;
  int squareBlockWidth  = 4;
  int squareBlockDepth  = 4;
  int numBlockLength = (squareLength / squareBlockLength);
  int numBlockWidth  = (squareWidth / squareBlockWidth);
  int numBlockDepth  = (squareDepth / squareBlockDepth);
  //numBlock = numBlockLength * numBlockWidth*numBlockDepth;
  //numBlock = 10003;
  int maxNumBlockVerts = 64;

  //FIMPtr->GraphPartition_Square(squareLength,squareWidth,squareDepth,
  //    squareBlockLength, squareBlockWidth, squareBlockDepth);    //use this for regular meshes
  FIMPtr->GraphPartition_METIS2( numBlock , maxNumBlockVerts); // use this for irregular meshes

  FIMPtr->m_numBlock = numBlock;
  FIMPtr->PartitionTets(numBlock);
  FIMPtr->GenerateData(maxIter,verbose);

  endtime = clock();
  double duration = (double)(endtime - starttime) *1000 / CLOCKS_PER_SEC;

  if (verbose)
    printf("Computing time : %.10lf ms\n",duration);

  return 0;
}

