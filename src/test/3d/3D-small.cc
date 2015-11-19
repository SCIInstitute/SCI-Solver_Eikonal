#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Small3D, EightTets) {
  Eikonal data(false);
  data.filename_ = TEST_DATA_DIR + std::string("sphere8");
  data.maxBlocks_ = 100;
  EXPECT_NO_THROW(data.solveEikonal());
  EXPECT_FLOAT_EQ(data.getFinalResult()[0],0.f);
  for (size_t i = 1; i < data.getFinalResult().size(); i ++) {
    EXPECT_FLOAT_EQ(data.getFinalResult()[i],19.58f);
  }
}
TEST(Small3D, ThirteenTets) {
  Eikonal data(false);
  data.filename_ = TEST_DATA_DIR + std::string("sphere13");
  data.maxBlocks_ = 100;
  EXPECT_NO_THROW(data.solveEikonal());
  EXPECT_NEAR(data.getFinalResult()[0],0.f,1e-4);
  for (size_t i = 1; i < data.getFinalResult().size(); i ++) {
    EXPECT_NEAR(data.getFinalResult()[i],19.58f,1e-4);
  }
}
