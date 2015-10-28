#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Converge3D, Unstructured) {
  Eikonal3D::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  EXPECT_NO_THROW(Eikonal3D::solveEikonal3D(data));
  for (size_t i = 2; i < Eikonal3D::numIterations(); i ++) {
    for (size_t j = 0; j < Eikonal3D::getResultAtIteration(i).size(); j ++) {
      float one = Eikonal3D::getResultAtIteration(i-1)[j];
      float zero = Eikonal3D::getResultAtIteration(i)[j];
      if (one == 0.) continue;
      EXPECT_TRUE(zero <= one);
    }
  }
  delete Eikonal3D::mesh_;
}
TEST(Converge3D, Structured) {
  Eikonal3D::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.isStructured_ = true;
  EXPECT_NO_THROW(Eikonal3D::solveEikonal3D(data));
  size_t sz = Eikonal3D::getFinalResult().size();
  for (size_t i = 2; i < Eikonal3D::numIterations(); i ++) {
    std::vector <float> resultA = Eikonal3D::getResultAtIteration(i-1);
    std::vector <float> resultB = Eikonal3D::getResultAtIteration(i);
    for (size_t j = 0; j < sz; j ++) {
      float one = resultA[j];
      float zero = resultB[j];
      if (one == 0.) continue;
      EXPECT_TRUE(zero <= one);
    }
  }
  delete Eikonal3D::mesh_;
}
