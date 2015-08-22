#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Basic3D, NonMaxUnstructured) {
  Eikonal::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  for (size_t i = 0; i < Eikonal::getFinalResult().size(); i ++) {
    EXPECT_TRUE(Eikonal::getFinalResult()[i] < 100.);
  }
  delete Eikonal::mesh_;
}
TEST(Basic3D, NonMaxStructured) {
  Eikonal::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.isStructured_ = true;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  std::vector<float> result = Eikonal::getFinalResult();
  for (size_t i = 0; i < result.size(); i ++) {
    EXPECT_TRUE(result[i] < 100.);
  }
  delete Eikonal::mesh_;
}
