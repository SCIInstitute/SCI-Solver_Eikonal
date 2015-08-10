#include "gtest/gtest.h"
#include "Eikonal2D.h"
TEST(Basic2D, NonMaxUnstructured) {
  Eikonal::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 1000;
  std::vector< std::vector <float> > results;
  EXPECT_NO_THROW((results = Eikonal::solveEikonal2D(data)));
  for (size_t i = 1; i < results.size(); i ++) {
    for (size_t j = 0; j < results[i].size(); j ++) {
      EXPECT_TRUE(results[i][j] <= results[i-1][j]);
    }
  }
}
TEST(Basic2D, NonMaxStructured) {
  Eikonal::Eikonal2D data;
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  std::vector< std::vector <float> > results;
  EXPECT_NO_THROW((results = Eikonal::solveEikonal2D(data)));
  for (size_t i = 1; i < results.size(); i ++) {
    for (size_t j = 0; j < results[i].size(); j ++) {
      EXPECT_TRUE(results[i][j] <= results[i-1][j]);
    }
  }
}
