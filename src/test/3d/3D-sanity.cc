#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Sanity3D, Filename) {
  Eikonal data(false);
  data.filename_ = "notarealfile";

  EXPECT_EXIT(data.solveEikonal(),
      ::testing::ExitedWithCode(1), "");
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  std::cout << data.filename_ << std::endl;
  EXPECT_NO_THROW(data.solveEikonal());
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  std::cout << data.filename_ << std::endl;
  data.isStructured_ = true;
  EXPECT_NO_THROW(data.solveEikonal());
}
