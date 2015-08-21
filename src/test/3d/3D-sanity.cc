#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Sanity3D, Filename) {
  Eikonal::Eikonal3D data;
  data.filename_ = "notarealfile";

  EXPECT_EXIT(Eikonal::solveEikonal3D(data),
      ::testing::ExitedWithCode(1), "");
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  delete Eikonal::mesh_;
}
