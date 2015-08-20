#ifndef __EIKONAL2D_H__
#define __EIKONAL2D_H__

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "meshFIM3d.h"
#include <math.h>
#include "tetgen.h"
#include "tetmesh.h"

#include <time.h>

namespace Eikonal {
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
        stopDistance_(50000.f),
        isStructured_(false),
        speedType_(ONE),
        squareLength_(16),
        squareWidth_(16),
        squareDepth_(16),
        squareBlockLength_(4),
        squareBlockWidth_(4),
        squareBlockDepth_(4),
        maxIterations_(1000)
    {}
      //3D data
      bool verbose_;
      std::string filename_;
      std::vector<int> seedPointList_;
      int maxBlocks_;
      int maxVertsPerBlock_;
      float stopDistance_;
      bool isStructured_;
      int speedType_; // ONE (1), CURVATURE (2), NOISE (3)
      int squareLength_, squareWidth_, squareDepth_;
      int squareBlockLength_, squareBlockWidth_, squareBlockDepth_;
      int maxIterations_;
  };

  //The static pointer to the mesh
  static TetMesh * mesh_ = NULL;

  /**
   * Creates the mesh, partitions the mesh, and runs the algorithm.
   *
   * @data The set of options for the Eikonal algorithm.
   *       The defaults are used if nothing is provided.
   * @return An array of iterations where each iteration has an array
   *         of values that correspond to vertex values at that iteration.
   */
  std::vector< std::vector< float> >
    solveEikonal3D(Eikonal3D data = Eikonal3D()) {

      clock_t starttime, endtime;
      starttime = clock ();
      tetgenio in;
      in.load_tetmesh((char*)data.filename_.c_str(),data.verbose_);

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

      if (!mesh_) exit(1);
      meshFIM3d FIMPtr;
      FIMPtr.SetSeedPoint(data.seedPointList_);
      FIMPtr.SetMesh(mesh_);//, data.speedType_);
      //FIMPtr.SetStopDistance(data.stopDistance_);//TODO
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

      //std::vector< std::vector< float > > results =
      FIMPtr.GenerateData(data.maxIterations_, data.verbose_);
      endtime = clock();
      double duration = (double)(endtime - starttime) * 1000/ CLOCKS_PER_SEC;

      if (data.verbose_)
        printf("Computing time : %.10lf ms\n",duration);

      //return results;
      return std::vector < std::vector < float > > ();
    }
}

#endif
