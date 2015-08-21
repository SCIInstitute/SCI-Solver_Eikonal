#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Converge3D, Unstructured) {
  Eikonal::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  for (size_t i = 2; i < Eikonal::numIterations(); i ++) {
    for (size_t j = 0; j < Eikonal::getResultAtIteration(i).size(); j ++) {
      float one = Eikonal::getResultAtIteration(i-1)[j];
      float zero = Eikonal::getResultAtIteration(i)[j];
      if (one == 0.) continue;
      EXPECT_TRUE(zero <= one);
    }
  }
  delete Eikonal::mesh_;
}
TEST(Converge3D, Structured) {
  Eikonal::Eikonal3D data;
  data.filename_ = TEST_DATA_DIR + std::string("CubeMesh_size16");
  data.squareLength_ = 16;
  data.squareWidth_ = 16;
  data.squareBlockLength_ = 4;
  data.squareBlockWidth_ = 4;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  for (size_t i = 2; i < Eikonal::numIterations(); i ++) {
    for (size_t j = 0; j < Eikonal::getResultAtIteration(i).size(); j ++) {
      float one = Eikonal::getResultAtIteration(i-1)[j];
      float zero = Eikonal::getResultAtIteration(i)[j];
      if (one == 0.) continue;
      EXPECT_TRUE(zero <= one);
    }
  }
  delete Eikonal::mesh_;
}
