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
      Eikonal2D(std::string fname = "../example_data/sphere_602verts.ply",
          bool verbose = false) :
        verbose_(verbose),
        filename_(fname),
        seedPointList_(std::vector<int>(1, 0)),
        maxBlocks_(10003),
        maxVertsPerBlock_(64),
        stopDistance_(50000.f),
        isStructured_(false),
        speedType_(ONE),
        squareLength_(1024),
        squareWidth_(1024),
        squareBlockLength_(8),
        squareBlockWidth_(8)
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
  };

  //The static pointer to the mesh 
  static TriMesh * mesh_ = NULL;
  
  /**
   * Creates the mesh, partitions the mesh, and runs the algorithm.
   *
   * @data The set of options for the Eikonal algorithm.
   *       The defaults are used if nothing is provided.
   * @return An array of iterations where each iteration has an array 
   *         of values that correspond to vertex values at that iteration.
   */
  std::vector< std::vector< float> >
    solveEikonal2D(Eikonal2D data = Eikonal2D()) {
      clock_t starttime, endtime;
      starttime = clock ();
      mesh_ = TriMesh::read(data.filename_.c_str(), data.verbose_);
      meshFIM2d FIMPtr;
      FIMPtr.SetSeedPoint(data.seedPointList_);
      FIMPtr.SetMesh(mesh_, data.speedType_);
      FIMPtr.SetStopDistance(data.stopDistance_);
      if (data.isStructured_)
        FIMPtr.GraphPartition_Square(data.squareLength_,data.squareWidth_,
            data.squareBlockLength_,
            data.squareBlockWidth_,
            data.verbose_);
      else
        FIMPtr.GraphPartition_METIS2( data.maxBlocks_,
            data.maxVertsPerBlock_, data.verbose_);

      FIMPtr.PartitionFaces(data.maxBlocks_);
      FIMPtr.InitializeLabels(data.maxBlocks_);

      std::vector< std::vector< float > > results =
        FIMPtr.GenerateData(data.maxBlocks_, data.verbose_);
      endtime = clock();
      double duration = (double)(endtime - starttime) * 1000/ CLOCKS_PER_SEC;

      if (data.verbose_)
        printf("Computing time : %.10lf ms\n",duration);

      return results;
    }
}

#endif
