#include "gtest/gtest.h"
#include "Eikonal2D.h"
TEST(Basic2D, NonMaxUnstructured) {
  Eikonal::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 11;
  data.stopDistance_ = 100.f;
  std::vector< std::vector <float> > results;
  EXPECT_NO_THROW((results = Eikonal::solveEikonal2D(data)));
  for (size_t i = 0; i < results[0].size(); i ++) {
    EXPECT_TRUE(results[results.size() - 1][i] < data.stopDistance_);
  }
}
TEST(Basic2D, NonMaxStructured) {
  Eikonal::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  data.stopDistance_ = 100.f;
  std::vector< std::vector <float> > results;
  EXPECT_NO_THROW((results = Eikonal::solveEikonal2D(data)));
  for (size_t i = 0; i < results[0].size(); i ++) {
    EXPECT_TRUE(results[results.size() - 1][i] < data.stopDistance_);
  }
}
