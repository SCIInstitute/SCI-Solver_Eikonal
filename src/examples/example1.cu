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
    if (strcmp(argv[i],"-v") == 0) {
      data.verbose_ = true;
    } else if (strcmp(argv[i],"-b") == 0) {
      if (i+1 >= argc) break;
      data.maxBlocks_ = atoi(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-i") == 0) {
      if (i+1 >= argc) break;
      data.filename_ = std::string(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-h") == 0) {
      printf("Usage: ./Example2 [OPTIONS]\n");
      printf("  -h            Show this help.\n");
      printf("  -v            Verbose output.\n");
      printf("  -i INPUT      Use this triangle mesh \n");
      printf("  -b MAX_BLOCKS Max # of blocks to use\n");
      exit(0);
    }
  data.solveEikonal();
  //write out the VTK files
  data.writeVTK();
  //we know that the solution should be the euclidean distance from the center.
  std::vector <float> solution;
  for (size_t i = 0; i < data.tetMesh_->vertices.size(); i++) {
    float x = data.tetMesh_->vertices[i][0];
    float y = data.tetMesh_->vertices[i][1];
    float z = data.tetMesh_->vertices[i][2];
    solution.push_back(std::sqrt((54.f - x)*(54.f-x)+(54.f-y)*(54.f-y)+(54.f-z)*(54.f-z)));
  }
  if (data.verbose_)
    data.printErrorGraph(solution);
  return 0;
}

