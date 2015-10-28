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
#include <Eikonal2D.h>

int main(int argc, char* argv[]) {
  Eikonal2D::Eikonal2D data;
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i],"-v") == 0) {
      data.verbose_ = true;
    } else if (strcmp(argv[i],"-i") == 0) {
      if (i+1 >= argc) break;
      data.filename_ = std::string(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-b") == 0) {
      if (i+1 >= argc) break;
      data.maxBlocks_ = atoi(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-h") == 0) {
      printf("Usage: ./Example2 [OPTIONS]\n");
      printf("  -h            Show this help.\n");
      printf("  -v            Verbose output.\n");
      printf("  -i INPUT      Use this triangle mesh \n");
      //# of blocks affects partitioning & convergence.
      //Adjust accordingly.
      printf("  -b MAX_BLOCKS Max # of blocks to use\n");
      exit(0);
    }
  Eikonal2D::solveEikonal2D(data);
  //write the output to file
  Eikonal2D::writeVTK();
  //the solution for the sphere examples (center 54,54,54, & radius 19.58)
  std::vector< float > solution;
  solution.resize(Eikonal2D::mesh_->vertices.size());
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = Eikonal2D::mesh_->vertices[i][0] - 54.f;
    float yDot = Eikonal2D::mesh_->vertices[i][1] - 54.f;
    float zDot = Eikonal2D::mesh_->vertices[i][2] - 54.f;
    solution[i] = 19.58f * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  if (data.verbose_)
    Eikonal2D::printErrorGraph(solution);
  return 0;
}
