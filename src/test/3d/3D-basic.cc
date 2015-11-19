#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Basic3D, NonMaxUnstructured) {
  Eikonal data(false);
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  EXPECT_NO_THROW(data.solveEikonal());
  for (size_t i = 0; i < data.getFinalResult().size(); i ++) {
    EXPECT_TRUE(data.getFinalResult()[i] < 100.);
  }
}
TEST(Basic3D, NonMaxStructured) {
  Eikonal data(false);
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.isStructured_ = true;
  EXPECT_NO_THROW(data.solveEikonal());
  std::vector<float> result = data.getFinalResult();
  for (size_t i = 0; i < result.size(); i ++) {
    EXPECT_TRUE(result[i] < 100.);
  }
}
