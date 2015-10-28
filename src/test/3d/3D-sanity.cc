#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Sanity3D, Filename) {
  Eikonal3D::Eikonal3D data;
  data.filename_ = "notarealfile";

  EXPECT_EXIT(Eikonal3D::solveEikonal3D(data),
      ::testing::ExitedWithCode(1), "");
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  EXPECT_NO_THROW(Eikonal3D::solveEikonal3D(data));
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.isStructured_ = true;
  EXPECT_NO_THROW(Eikonal3D::solveEikonal3D(data));
  delete Eikonal3D::mesh_;
}
