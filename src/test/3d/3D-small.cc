#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Small3D, EightTets) {
  Eikonal3D::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere8");
  data.maxBlocks_ = 100;
  EXPECT_NO_THROW(Eikonal3D::solveEikonal3D(data));
  EXPECT_FLOAT_EQ(Eikonal3D::getFinalResult()[0],0.f);
  for (size_t i = 1; i < Eikonal3D::getFinalResult().size(); i ++) {
    EXPECT_FLOAT_EQ(Eikonal3D::getFinalResult()[i],19.58f);
  }
  delete Eikonal3D::mesh_;
}
TEST(Small3D, ThirteenTets) {
  Eikonal3D::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere13");
  data.maxBlocks_ = 100;
  EXPECT_NO_THROW(Eikonal3D::solveEikonal3D(data));
  EXPECT_NEAR(Eikonal3D::getFinalResult()[0],0.f,1e-4);
  for (size_t i = 1; i < Eikonal3D::getFinalResult().size(); i ++) {
    EXPECT_NEAR(Eikonal3D::getFinalResult()[i],19.58f,1e-4);
  }
  delete Eikonal3D::mesh_;
}
