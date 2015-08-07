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
#include <Eikonal2D.h>

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
      exit(0);
    }

  Eikonal::Eikonal2D data(filename);
  data.verbose_ = verbose;

  std::vector< std::vector <float> >
    results = Eikonal::solveEikonal2D(data);

  // find the analytical solution to each vertex and compare.
  std::vector< float > solution;
  solution.resize(Eikonal::mesh_->vertices.size());
  float radius = 19.58f; //we know the radius of these spheres.
  std::vector<float> center;
  //we know the center of these spheres.
  center.push_back(54.f);
  center.push_back(54.f);
  center.push_back(54.f);
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = Eikonal::mesh_->vertices[i][0] - center[0];
    float yDot = Eikonal::mesh_->vertices[i][1] - center[1];
    float zDot = Eikonal::mesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for each iteration
  std::vector<float> rmsError;
  rmsError.resize(results.size());
  for (size_t i = 0; i < results.size(); i++) {
    float sum = 0.f;
    for (size_t j = 0; j < solution.size(); j++) {
      float err = std::abs(solution[j] - results[i][j]);
      sum +=  err * err;
    }
    rmsError[i] = std::sqrt(sum / static_cast<float>(solution.size()));
  }
  //determine the log range
  float max_err = rmsError[0];
  float min_err = rmsError[rmsError.size() - 1];
  int max_log = -10, min_log = 10;
  while (std::pow(10,max_log) < max_err) max_log++;
  while (std::pow(10,min_log) > min_err) min_log--;
  // print the error graph
  printf("\n\nlog(Err)|\n");
  bool printTick = true;
  for(int i = max_log ; i >= min_log; i--) {
    if (printTick) {
      printf("   10^%2d|",i);
    } else {
      printf("        |");
    }
    for (size_t j = 0; j < results.size(); j++) {
      if (rmsError[j] > std::pow(10,i) &&
          rmsError[j] < std::pow(10,i+1))
        printf("*");
      else
        printf(" ");
    }
    printf("\n");
    printTick = !printTick;
  }
  printf("--------|------------------------------------------");
  printf("  Converged to: %.4f\n",rmsError[rmsError.size() - 1]);
  printf("        |0    5    10    15    20    25    30    35\n");
  printf("                   Iteration\n");

  return 0;
}
