#ifndef __EIKONAL3D_H__
#define __EIKONAL3D_H__

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "meshFIM3d.h"
#include <math.h>
#include "tetgen.h"
#include "tetmesh.h"

#include <time.h>

namespace Eikonal3D {
  /** The class that represents all of the available options for Eikonal3D */
  class Eikonal3D {
    public:
      Eikonal3D(std::string fname = "../src/test/test_data/sphere8",
          bool verbose = false) :
        verbose_(verbose),
        filename_(fname),
        seedPointList_(std::vector<int>(1, 0)),
        maxBlocks_(1000),
        maxVertsPerBlock_(64),
        isStructured_(false),
        squareLength_(16),
        squareWidth_(16),
        squareDepth_(16),
        squareBlockLength_(1),
        squareBlockWidth_(1),
        squareBlockDepth_(1),
        maxIterations_(1000)

    {}
      //3D data
      bool verbose_;
      std::string filename_;
      std::vector<int> seedPointList_;
      int maxBlocks_;
      int maxVertsPerBlock_;
      bool isStructured_;
      int squareLength_, squareWidth_, squareDepth_;
      int squareBlockLength_, squareBlockWidth_, squareBlockDepth_;
      int maxIterations_;
  };

  //The static pointer to the mesh
  static TetMesh * mesh_ = NULL;
  //the answer vector
  static std::vector < std::vector <float> > iteration_values_;
  //accessor functions to the results.
  std::vector < float >& getFinalResult() {
    return iteration_values_.at(iteration_values_.size() - 1);
  }
  std::vector < float >& getResultAtIteration(size_t i) {
    return iteration_values_.at(i);
  }
  size_t numIterations() { return iteration_values_.size(); }
  void writeVTK() {
    meshFIM3d FIMPtr;
    FIMPtr.SetMesh(mesh_);
    FIMPtr.writeVTK(iteration_values_);
  }

  /**
   * Creates the mesh, partitions the mesh, and runs the algorithm.
   *
   * @data The set of options for the Eikonal algorithm.
   *       The defaults are used if nothing is provided.
   */
  void solveEikonal3D(Eikonal3D data = Eikonal3D()) {

    clock_t starttime, endtime;
    starttime = clock ();
    tetgenio in;
    if (!in.load_tetmesh((char*)data.filename_.c_str(),data.verbose_))
      exit(1);

    mesh_ = new TetMesh();
    mesh_->init(
        in.pointlist,
        in.numberofpoints,
        in.trifacelist,
        in.numberoffacets,
        in.tetrahedronlist,
        in.numberoftetrahedra,
        in.numberoftetrahedronattributes,
        in.tetrahedronattributelist, data.verbose_ );
    mesh_->need_neighbors(data.verbose_);
    mesh_->need_adjacenttets(data.verbose_);
    mesh_->need_tet_virtual_tets(data.verbose_);

    meshFIM3d FIMPtr;
    FIMPtr.SetSeedPoint(data.seedPointList_);
    FIMPtr.SetMesh(mesh_);
    FIMPtr.InitSpeedMat();
    if (data.isStructured_) {
      int numBlockLength = (data.squareLength_ / data.squareBlockLength_);
      int numBlockWidth  = (data.squareWidth_ / data.squareBlockWidth_);
      int numBlockDepth  = (data.squareDepth_ / data.squareBlockDepth_);
      data.maxBlocks_ = numBlockLength * numBlockWidth * numBlockDepth;
      FIMPtr.GraphPartition_Square(
          data.squareLength_, data.squareWidth_,data.squareDepth_,
          data.squareBlockLength_,
          data.squareBlockWidth_, data.squareBlockDepth_,
          data.verbose_);

    } else {
      FIMPtr.GraphPartition_METIS2( data.maxBlocks_,
          data.maxVertsPerBlock_, data.verbose_);
    }
    FIMPtr.m_numBlock = data.maxBlocks_;
    FIMPtr.PartitionTets(data.maxBlocks_,data.verbose_);
    iteration_values_ =
      FIMPtr.GenerateData(data.maxIterations_, data.verbose_);
    endtime = clock();
    double duration = (double)(endtime - starttime) * 1000/ CLOCKS_PER_SEC;

    if (data.verbose_)
      printf("Computing time : %.10lf ms\n",duration);
  }

  /**
   * This function uses the provided analytical solutions to
   * visualize the algorithm's error after each iteration.
   *
   * @param solution The vector of expected solutions.
   */
  void printErrorGraph(std::vector<float> solution) {

    // now calculate the RMS error for each iteration
    std::vector<float> rmsError;
    rmsError.resize(numIterations());
    for (size_t i = 0; i < numIterations(); i++) {
      float sum = 0.f;
      std::vector<float> result = getResultAtIteration(i);
      for (size_t j = 0; j < solution.size(); j++) {
        float err = std::abs(solution[j] - result[j]);
        sum +=  err * err;
      }
      rmsError[i] = std::sqrt(sum / static_cast<float>(solution.size()));
    }
    //determine the log range
    float max_err = rmsError[0];
    float min_err = rmsError[rmsError.size() - 1];
    int max_log = -10, min_log = 10;
    while (std::pow(static_cast<float>(10),max_log) < max_err) max_log++;
    while (std::pow(static_cast<float>(10),min_log) > min_err) min_log--;
    // print the error graph
    printf("\n\nlog(Err)|\n");
    bool printTick = true;
    for(int i = max_log ; i >= min_log; i--) {
      if (printTick) {
        printf("   10^%2d|",i);
      } else {
        printf("        |");
      }
      for (size_t j = 0; j < numIterations(); j++) {
        if (rmsError[j] > std::pow(static_cast<float>(10),i) &&
            rmsError[j] < std::pow(static_cast<float>(10),i+1))
          printf("*");
        else
          printf(" ");
      }
      printf("\n");
      printTick = !printTick;
    }
    printf("--------|------------------------------------------");
    printf("  Converged to: %.4f\n",rmsError[rmsError.size() - 1]);
    printf("        |1   5    10   15   20   25   30   35\n");
    printf("                   Iteration\n");
  }
}

#endif
