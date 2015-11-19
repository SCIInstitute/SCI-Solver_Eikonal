#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Sanity2D, Filename) {
  Eikonal data(true);
  data.filename_ = "notarealfile";

  EXPECT_EXIT(data.solveEikonal(),
      ::testing::ExitedWithCode(1), "");
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(data.solveEikonal());
  data.filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(data.solveEikonal());
}
