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

int main(int argc, char *argv[])
{
  //3D options
  Eikonal data(false);
  //input filename (minus extension)
  data.filename_ = "../src/test/test_data/sphere334";
  for (int i = 0; i < argc; i++)
    if (strcmp(argv[i], "-v") == 0) {
    data.verbose_ = true;
    } else if (strcmp(argv[i], "-m") == 0) {
    if (i + 1 >= argc) break;
    data.maxIterations_ = atoi(argv[i + 1]);
    i++;
    } else if (strcmp(argv[i], "-b") == 0) {
      if (i + 1 >= argc) break;
      data.maxBlocks_ = atoi(argv[i + 1]);
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
    } else if (strcmp(argv[i], "-i") == 0) {
      if (i + 1 >= argc) break;
      data.filename_ = std::string(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-h") == 0) {
      printf("Usage: ./Example1 [OPTIONS]\n");
      printf("  -h              Show this help.\n");
      printf("  -v              Verbose output.\n");
      printf("  -i INPUT        Use this tet mesh \n");
      printf("  -b MAX_BLOCKS   Max # of blocks to use\n");
      printf("  -m MAX_ITER     Max # of iterations before quit\n");
      // The tensors are 6 unique values in a matrix per tet
      // [ 0 1 2 ]
      // [ 1 3 4 ]
      // [ 2 4 5 ]
      printf("  -x MATRIX_FILE  File of tensor matrices per tet (N*6 floats).\n");
      printf("                  N floats if considered scalar speeds per tet.\n");
      exit(0);
    }
    data.solveEikonal();
    //write out the VTK files
    data.writeVTK(false); //true to output values at each iter.
    //we know that the solution should be the euclidean distance from the center.
    std::vector <float> solution;
    for (size_t i = 0; i < data.tetMesh_->vertices.size(); i++) {
      float x = data.tetMesh_->vertices[i][0];
      float y = data.tetMesh_->vertices[i][1];
      float z = data.tetMesh_->vertices[i][2];
      solution.push_back(std::sqrt((0.f - x)*(0.f - x) + (0.f - y)*(0.f - y) + (0.f - z)*(0.f - z)));
    }
    if (data.verbose_)
      data.printErrorGraph(solution);
    return 0;
}

