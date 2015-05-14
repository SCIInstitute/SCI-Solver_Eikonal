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
#include "TriMesh.h"
#include <cstring>
#include "meshFIM2d.h"
#include <time.h>

/************************************************************************/
/* main                                                           */
/************************************************************************/
int main(int argc, char* argv[]) {
  //Verbose option
  bool verbose = false;
  //input filename (minus extension)
  std::string filename;
  for (int i = 0; i < argc; i++)
    if (strcmp(argv[i],"-v") == 0) {
      verbose = true;
    } else if (strcmp(argv[i],"-i") == 0) {
      if (i+1 >= argc) break;
      filename = std::string(argv[i+1]);
      i++;
    }
  if (filename.empty())
    filename = "../example_data/SquareMesh_size16.ply";
  clock_t starttime, endtime;

  TriMesh *themesh = TriMesh::read(filename.c_str());

  //themesh->need_normals();
  //themesh->need_tstrips();
  themesh->need_bsphere();
  //themesh->need_faceedges();
  //themesh->need_across_edge();

  meshFIM2d* FIMPtr = new meshFIM2d;

  starttime = clock ();

  std::vector<int> seedPointList(1,0/*,currentVert*/);

  //int squareLength = 1024;
  //int squareWidth = 1024;
  //int squareBlockLength = 8;
  //int squareBlockWidth  = 8;
  //int numBlockLength = (squareLength / squareBlockLength);
  //int numBlockWidth  = (squareWidth / squareBlockWidth);
  //int numBlock = numBlockLength * numBlockWidth;

  int numBlock = /*0*//*16226*//*18729*//*16624*//*5210*//*2346*//*803*//*14295*//*1600*/10003 /*2557*//*3487*//*1566*//*177*//*950*//*3487*//*1750*//*1150*//*2309*/ /*175*//*2800*/;  //for 59021verts, 950 for 64; for 72202verts, 1150 for 64, 2309 for 32.; for dragon, 1600 for 64,for dragon iso, 3487 for 32; for 98687, 1566 for 64; for dragon.ts, 2557 for 64; for dragon.ts_maxSF0.5, 10003 for 64; for square.1.ply, 14295 for 64;for sphereR40_iso.ply, 803 for 64; 2346 for sphereR60_147237.ply;5210 for square328k; 16624 for square_size1024;18729 for square_1.1m;16226 for sphereR60_1024k_split and 1024k_256split
  int maxNumBlockVerts = 64;

  FIMPtr->SetSeedPoint(seedPointList);
  FIMPtr->SetMesh(themesh);
  FIMPtr->SetStopDistance(50000.0);
  FIMPtr->GraphPartition_METIS2( numBlock, maxNumBlockVerts);

  FIMPtr->PartitionFaces(numBlock);
  FIMPtr->InitializeLabels(numBlock);

  FIMPtr->GenerateData(numBlock);

  endtime = clock();
  double duration = (double)(endtime - starttime) * 1000/ CLOCKS_PER_SEC;

  if (verbose)
    printf("Computing time : %.10lf ms\n",duration);

  delete FIMPtr;
  return 0;
}
