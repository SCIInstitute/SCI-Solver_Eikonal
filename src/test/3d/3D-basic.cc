#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Basic3D, NonMaxUnstructured) {
  Eikonal::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  data.maxBlocks_ = 100;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  for (size_t i = 0; i < Eikonal::getFinalResult().size(); i ++) {
    EXPECT_TRUE(Eikonal::getFinalResult()[i] < 100.);
  }
  delete Eikonal::mesh_;
}
TEST(Basic3D, NonMaxStructured) {
  Eikonal::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  for (size_t i = 0; i < Eikonal::getFinalResult().size(); i ++) {
    EXPECT_TRUE(Eikonal::getFinalResult()[i] < 100.);
  }
  delete Eikonal::mesh_;
}
