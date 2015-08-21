#include "gtest/gtest.h"
#include "Eikonal2D.h"
TEST(Converge2D, Unstructured) {
  Eikonal::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(Eikonal::solveEikonal2D(data));
  for (size_t i = 1; i < Eikonal::numIterations(); i ++) {
    for (size_t j = 0; j < Eikonal::getFinalResult().size(); j ++) {
      EXPECT_TRUE(Eikonal::getResultAtIteration(i)[j] <= 
          Eikonal::getResultAtIteration(i)[j]);
    }
  }
}
TEST(Converge2D, Structured) {
  Eikonal::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(Eikonal::solveEikonal2D(data));
  for (size_t i = 1; i < Eikonal::numIterations(); i ++) {
    for (size_t j = 0; j < Eikonal::getFinalResult().size(); j ++) {
      EXPECT_TRUE(Eikonal::getResultAtIteration(i)[j] <= 
          Eikonal::getResultAtIteration(i)[j]);
    }
  }
}
