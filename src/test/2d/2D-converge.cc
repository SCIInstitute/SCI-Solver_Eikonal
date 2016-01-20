#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Converge2D, Unstructured) {
  Eikonal data(true);
  data.filename_ = TEST_DATA_DIR + std::string("sphere_74verts.ply");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(data.solveEikonal());
  for (size_t i = 1; i < data.numIterations(); i ++) {
    for (size_t j = 0; j < data.getFinalResult().size(); j ++) {
      EXPECT_TRUE(data.getResultAtIteration(i)[j] <= 
          data.getResultAtIteration(i)[j]);
    }
  }
}
TEST(Converge2D, Structured) {
  Eikonal data(true);
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(data.solveEikonal());
  for (size_t i = 1; i < data.numIterations(); i ++) {
    for (size_t j = 0; j < data.getFinalResult().size(); j ++) {
      EXPECT_TRUE(data.getResultAtIteration(i)[j] <= 
          data.getResultAtIteration(i)[j]);
    }
  }
}
