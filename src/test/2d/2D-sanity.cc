#include "gtest/gtest.h"
#include "Eikonal2D.h"
TEST(Sanity2D, Filename) {
  Eikonal2D::Eikonal2D data;
  data.filename_ = "notarealfile";

  EXPECT_EXIT(Eikonal2D::solveEikonal2D(data),
      ::testing::ExitedWithCode(1), "");
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
}
