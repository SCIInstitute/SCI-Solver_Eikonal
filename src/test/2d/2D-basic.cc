#include "gtest/gtest.h"
#include "Eikonal2D.h"
TEST(Basic2D, NonMaxUnstructured) {
  Eikonal2D::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 11;
  data.stopDistance_ = 100.f;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  for (size_t i = 0; i < Eikonal2D::getFinalResult().size(); i ++) {
    EXPECT_TRUE(Eikonal2D::getFinalResult()[i] < data.stopDistance_);
  }
}
TEST(Basic2D, NonMaxStructured) {
  Eikonal2D::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  data.stopDistance_ = 100.f;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  for (size_t i = 0; i < Eikonal2D::getFinalResult().size(); i ++) {
    EXPECT_TRUE(Eikonal2D::getFinalResult()[i] < data.stopDistance_);
  }
}
