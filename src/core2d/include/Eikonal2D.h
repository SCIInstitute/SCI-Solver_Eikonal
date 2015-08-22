#ifndef __EIKONAL2D_H__
#define __EIKONAL2D_H__

#include "TriMesh.h"
#include "meshFIM2d.h"
#include <cstring>
#include <time.h>

namespace Eikonal {
  /** The class that represents all of the available options for Eikonal2D */
  class Eikonal2D {
    public:
      Eikonal2D(std::string fname = "../src/test/test_data/sphere_266verts.ply",
          bool verbose = false) :
        verbose_(verbose),
        filename_(fname),
        seedPointList_(std::vector<int>(1, 0)),
        maxBlocks_(10003),
        maxVertsPerBlock_(64),
        stopDistance_(50000.f),
        isStructured_(false),
        speedType_(ONE),
        squareLength_(16),
        squareWidth_(16),
        squareBlockLength_(1),
        squareBlockWidth_(1),
        maxIterations_(1000)
    {}
      //2D data
      bool verbose_;
      std::string filename_;
      std::vector<int> seedPointList_;
      int maxBlocks_;
      int maxVertsPerBlock_;
      float stopDistance_;
      bool isStructured_;
      int speedType_; // ONE (1), CURVATURE (2), NOISE (3)
      int squareLength_, squareWidth_;
      int squareBlockLength_, squareBlockWidth_;
      int maxIterations_;
  };

  //The static pointer to the mesh
  static TriMesh * mesh_ = NULL;
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

  /**
   * Creates the mesh, partitions the mesh, and runs the algorithm.
   *
   * @data The set of options for the Eikonal algorithm.
   *       The defaults are used if nothing is provided.
   */
  void solveEikonal2D(Eikonal2D data = Eikonal2D()) {
    clock_t starttime, endtime;
    starttime = clock ();
    mesh_ = TriMesh::read(data.filename_.c_str(), data.verbose_);
    if (!mesh_) exit(1);
    meshFIM2d FIMPtr;
    FIMPtr.SetSeedPoint(data.seedPointList_);
    FIMPtr.SetMesh(mesh_, data.speedType_);
    FIMPtr.SetStopDistance(data.stopDistance_);
    if (data.isStructured_) {
      int numBlockLength = (data.squareLength_ / data.squareBlockLength_);
      int numBlockWidth  = (data.squareWidth_ / data.squareBlockWidth_);
      data.maxBlocks_ = numBlockLength * numBlockWidth;
      FIMPtr.GraphPartition_Square(data.squareLength_,data.squareWidth_,
          data.squareBlockLength_,
          data.squareBlockWidth_,
          data.verbose_);

    } else {
      FIMPtr.GraphPartition_METIS2( data.maxBlocks_,
          data.maxVertsPerBlock_, data.verbose_);
    }
    FIMPtr.PartitionFaces(data.maxBlocks_);
    FIMPtr.InitializeLabels(data.maxBlocks_);

    iteration_values_ =
      FIMPtr.GenerateData(data.maxBlocks_, data.maxIterations_, data.verbose_);
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
