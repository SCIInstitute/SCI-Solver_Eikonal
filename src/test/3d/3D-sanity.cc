#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Sanity3D, Filename) {
  Eikonal data(false);
  data.filename_ = "notarealfile";

  EXPECT_EXIT(data.solveEikonal(),
      ::testing::ExitedWithCode(1), "");
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  EXPECT_NO_THROW(data.solveEikonal());
  delete data.tetMesh_;
  data.tetMesh_ = NULL;
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.isStructured_ = true;
  EXPECT_NO_THROW(data.solveEikonal());
}
