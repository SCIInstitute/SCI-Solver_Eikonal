#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Basic2D, NonMaxUnstructured) {
  Eikonal data(true);
  data.filename_ = TEST_DATA_DIR + std::string("sphere_74verts.ply");
  data.maxBlocks_ = 11;
  data.stopDistance_ = 700.f;
  EXPECT_NO_THROW(data.solveEikonal());
  for (size_t i = 0; i < data.getFinalResult().size(); i ++) {
    EXPECT_TRUE(data.getFinalResult()[i] < data.stopDistance_);
  }
}
TEST(Basic2D, NonMaxStructured) {
  Eikonal data(true);
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  data.stopDistance_ = 100.f;
  data.isStructured_ = true;
  EXPECT_NO_THROW(data.solveEikonal());
  for (size_t i = 0; i < data.getFinalResult().size(); i ++) {
    EXPECT_TRUE(data.getFinalResult()[i] < data.stopDistance_);
  }
}
