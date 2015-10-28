#include "gtest/gtest.h"
#include "Eikonal2D.h"
TEST(Converge2D, Unstructured) {
  Eikonal2D::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  for (size_t i = 1; i < Eikonal2D::numIterations(); i ++) {
    for (size_t j = 0; j < Eikonal2D::getFinalResult().size(); j ++) {
      EXPECT_TRUE(Eikonal2D::getResultAtIteration(i)[j] <= 
          Eikonal2D::getResultAtIteration(i)[j]);
    }
  }
}
TEST(Converge2D, Structured) {
  Eikonal2D::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  for (size_t i = 1; i < Eikonal2D::numIterations(); i ++) {
    for (size_t j = 0; j < Eikonal2D::getFinalResult().size(); j ++) {
      EXPECT_TRUE(Eikonal2D::getResultAtIteration(i)[j] <= 
          Eikonal2D::getResultAtIteration(i)[j]);
    }
  }
}
