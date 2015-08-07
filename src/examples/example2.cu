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
#include <cstring>
#include <time.h>
#include "TriMesh.h"
#include "meshFIM2d.h"

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
    } else if (strcmp(argv[i],"-h") == 0) {
      printf("Usage: ./Example2 [OPTIONS]\n");
      printf("  -h            Show this help.\n");
      printf("  -v            Verbose output.\n");
      printf("  -i INPUT      Use this triangle mesh \n");
      printf("                [default ../example_data/sphere_4196verts.ply]\n");
      exit(0);
    }
  if (filename.empty())
    filename = "../example_data/sphere_4196verts.ply";
  clock_t starttime, endtime;

  TriMesh *themesh = TriMesh::read(filename.c_str(), verbose);

  meshFIM2d* FIMPtr = new meshFIM2d;

  starttime = clock ();

  std::vector<int> seedPointList(1,0/*,currentVert*/);

  int numBlock = 10003;
  int maxNumBlockVerts = 64;
  float stopDist = 50000.f;

  FIMPtr->SetSeedPoint(seedPointList);
  FIMPtr->SetMesh(themesh);
  FIMPtr->SetStopDistance(stopDist);
  FIMPtr->GraphPartition_METIS2( numBlock, maxNumBlockVerts, verbose);

  FIMPtr->PartitionFaces(numBlock);
  FIMPtr->InitializeLabels(numBlock);

  std::vector< std::vector< float > > results =
    FIMPtr->GenerateData(numBlock, verbose);

  endtime = clock();
  double duration = (double)(endtime - starttime) * 1000/ CLOCKS_PER_SEC;

  if (verbose)
    printf("Computing time : %.10lf ms\n",duration);

  delete FIMPtr;
  // find the analytical solution to each vertex and compare.
  std::vector< float > solution;
  solution.resize(themesh->vertices.size());
  float radius = 19.58f; //we know the radius of these spheres.
  std::vector<float> center;
  //we know the center of these spheres.
  center.push_back(54.f);
  center.push_back(54.f);
  center.push_back(54.f);
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = themesh->vertices[i][0] - center[0];
    float yDot = themesh->vertices[i][1] - center[1];
    float zDot = themesh->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for each iteration
  std::vector<float> rmsError;
  rmsError.resize(results.size());
  size_t badVerts = 0;
  for (size_t i = 0; i < results.size(); i++) {
    float sum = 0.f;
    for (size_t j = 0; j < solution.size(); j++) {
      float err = std::abs(solution[j] - results[i][j]);
      if (i == results.size() - 1 && err > stopDist) {
        badVerts++;
        err = 0.f;
      }
      sum +=  err * err;
    }
    rmsError[i] = std::sqrt(sum / static_cast<float>(solution.size()));
    std::cout << rmsError[i] << std::endl;

  }
  std::cout << "# Bad vertex values: " << badVerts << std::endl;
  return 0;
}
