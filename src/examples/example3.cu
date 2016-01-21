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

#include <Eikonal.h>
//this example prepares for a "CURVATURE" eikonal solution.
int main(int argc, char *argv[])
{
  std::string fname, type;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {
      if (i + 1 >= argc) break;
      fname = std::string(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-s") == 0) {
      if (i + 1 >= argc) break;
      type = std::string(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-h") == 0) {
      printf("Usage: ./Example3 [OPTIONS]\n");
      printf("  -i INPUT      Use this triangle mesh \n");
      printf("  -s SPEEDTYPE    Speed type is [ONE], CURVATURE, or NOISE.\n");
      exit(0);
    }
  }
  bool triMesh = fname.rfind(".ply") != std::string::npos;
  Eikonal data(triMesh);
  data.verbose_ = true;
  data.filename_ = fname;
  if (type == "CURVATURE") {
    data.speedType_ = CURVATURE;
  } else if (type == "NOISE") {
    data.speedType_ = NOISE;
  } else {
    data.speedType_ = ONE;
  }
  data.initializeMesh();
  // initialize exactly 1 unit in from cube/square as seed points. (R3 {x,y,z} : 0 -> 15)
  std::vector<float> vals;
  size_t num = triMesh ? data.triMesh_->vertices.size() : data.tetMesh_->vertices.size();
  for (size_t i = 0; i < num; i++) {
    bool found = false;
    for (size_t j = 0; j < 3; j++) {
      if (triMesh) {
        found |= data.triMesh_->vertices[i][j] == 1.;
        found |= data.triMesh_->vertices[i][j] == 14.;
      } else {
        found |= data.tetMesh_->vertices[i][j] == 1.;
        found |= data.tetMesh_->vertices[i][j] == 14.;
      }
    }
    vals.push_back(found ? 0.f : LARGENUM);
  }
  data.initializeVertices(vals);
  data.solveEikonal();
  //write out the VTK files
  data.writeVTK(false); //true to output values at each iter.
  return 0;
}

