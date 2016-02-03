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

int main(int argc, char* argv[]) {
  Eikonal data(true);
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      data.verbose_ = true;
    } else if (strcmp(argv[i], "-i") == 0) {
      if (i + 1 >= argc) break;
      data.filename_ = std::string(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-s") == 0) {
      if (i + 1 >= argc) break;
      std::string type = std::string(argv[i + 1]);
      if (type == "CURVATURE") {
        data.speedType_ = CURVATURE;
      } else if (type == "NOISE") {
        data.speedType_ = NOISE;
      } else if (type == "MANUAL"){
        data.speedType_ = MANUAL_SPEED;
      }
      i++;
    } else if (strcmp(argv[i], "-x") == 0) {
      while (i + 1 < argc && argv[i + 1][0] != '-') {
        std::ifstream mat(argv[++i]);
        while (mat.good()) {
          float val;
          mat >> val;
          if (!mat.good()) break;
          data.speedMtxMultipliers_.push_back(val);
        }
        mat.close();
      }
    } else if (strcmp(argv[i], "-b") == 0) {
      if (i + 1 >= argc) break;
      data.maxVertsPerBlock_ = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-n") == 0) {
      if (i + 1 >= argc) break;
      data.maxIterations_= atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-h") == 0) {
      printf("Usage: ./Example2 [OPTIONS]\n");
      printf("  -h              Show this help.\n");
      printf("  -v              Verbose output.\n");
      printf("  -i INPUT        Use this triangle mesh \n");
      //# of verts/block affects partitioning & convergence.
      //Adjust accordingly.
      printf("  -b MAX_BLK_VERT Max # of verts/block to use\n");
      printf("  -n MAX_ITER     Max # of iterations to run\n");
      printf("  -s SPEEDTYPE    Speed type is [ONE], CURVATURE, or NOISE.\n");
      printf("  -x SCALAR_FILE  File with speed scalars per face.\n");
      exit(0);
    }
  }
  data.initializeMesh();
  std::ofstream out("square_scalars.txt");
  for (size_t i = 0; i < data.triMesh_->faces.size(); i++) {
    point p = (data.triMesh_->vertices[data.triMesh_->faces[i][0]] +
      data.triMesh_->vertices[data.triMesh_->faces[i][1]] +
      data.triMesh_->vertices[data.triMesh_->faces[i][2]]) / 3.f;
    if (len(p - point(16., 16., 0.)) < 6.) {
      out << 20. << std::endl;
    } else {
      out << 1. << std::endl;
    }
  }
  out.close(); 

  data.solveEikonal();
  //write the output to file
  data.writeVTK(false);
  //the solution for the sphere examples (center 0,0,0, & radius 100)
  std::vector< float > solution;
  solution.resize(data.triMesh_->vertices.size());
  point first = data.triMesh_->vertices[0];
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data.triMesh_->vertices[i][0];
    float yDot = data.triMesh_->vertices[i][1];
    float zDot = data.triMesh_->vertices[i][2];
    solution[i] = 100.f * std::acos(
      (first[0] * xDot + first[1] * yDot + first[2] * zDot) /
      std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot) /
      len(first));
  }
  if (data.verbose_)
    data.printErrorGraph(solution);
  return 0;
}
