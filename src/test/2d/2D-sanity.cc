#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Sanity2D, Filename) {
  Eikonal* data = new Eikonal(true);
  data->filename_ = "notarealfile";

  EXPECT_EXIT(data->solveEikonal(),
      ::testing::ExitedWithCode(1), "");
  delete data;
  data = new Eikonal(true);
  data->filename_ = TEST_DATA_DIR + std::string("sphere_74verts.ply");
  EXPECT_NO_THROW(data->solveEikonal());
  delete data;
  data = new Eikonal(true);
  data->filename_ = TEST_DATA_DIR + std::string("SquareMesh_size16.ply");
  data->squareLength_ = 16;
  data->squareWidth_ = 16;
  data->squareBlockLength_ = 4;
  data->squareBlockWidth_ = 4;
  data->isStructured_ = true;
  EXPECT_NO_THROW(data->solveEikonal());
  delete data;
}
