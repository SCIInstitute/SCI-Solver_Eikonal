#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Converge3D, Unstructured) {
  Eikonal data(false);
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  EXPECT_NO_THROW(data.solveEikonal());
  for (size_t i = 2; i < data.numIterations(); i ++) {
    for (size_t j = 0; j < data.getResultAtIteration(i).size(); j ++) {
      float one = data.getResultAtIteration(i-1)[j];
      float zero = data.getResultAtIteration(i)[j];
      if (one == 0.) continue;
      EXPECT_TRUE(zero <= one);
    }
  }
}
TEST(Converge3D, Structured) {
  Eikonal data(false);
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.isStructured_ = true;
  EXPECT_NO_THROW(data.solveEikonal());
  size_t sz = data.getFinalResult().size();
  for (size_t i = 2; i < data.numIterations(); i ++) {
    std::vector <float> resultA = data.getResultAtIteration(i-1);
    std::vector <float> resultB = data.getResultAtIteration(i);
    for (size_t j = 0; j < sz; j ++) {
      float one = resultA[j];
      float zero = resultB[j];
      if (one == 0.) continue;
      EXPECT_TRUE(zero <= one);
    }
  }
}
