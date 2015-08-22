#include "gtest/gtest.h"
#include "Eikonal3D.h"
TEST(Improve3D, DecreaseRMSError) {
  Eikonal::Eikonal3D data;
  std::vector<float> rmsError;
  //low granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere334");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  // find the analytical solution to each vertex and compare.
  std::vector< float > solution;
  solution.resize(Eikonal::mesh_->vertices.size());
  for(size_t i = 0; i < Eikonal::mesh_->vertices.size(); i++) {
    float x = Eikonal::mesh_->vertices[i][0];
    float y = Eikonal::mesh_->vertices[i][1];
    float z = Eikonal::mesh_->vertices[i][2];
    solution[i] = std::sqrt((54. - x)*(54.-x)+(54.-y)*(54.-y)+(54.-z)*(54.-z));
  }
  // now calculate the RMS error for this run
  float sum = 0.f;
  std::vector <float> result = Eikonal::getFinalResult();
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - result[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  delete Eikonal::mesh_;
  //mid granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere8092");
  data.maxBlocks_ = 1000;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  // find the analytical solution to each vertex and compare.
  solution.resize(Eikonal::mesh_->vertices.size());
  for(size_t i = 0; i < Eikonal::mesh_->vertices.size(); i++) {
    float x = Eikonal::mesh_->vertices[i][0];
    float y = Eikonal::mesh_->vertices[i][1];
    float z = Eikonal::mesh_->vertices[i][2];
    solution[i] = std::sqrt((54. - x)*(54.-x)+(54.-y)*(54.-y)+(54.-z)*(54.-z));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  result = Eikonal::getFinalResult();
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - result[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  delete Eikonal::mesh_;
  //highest granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere19499");
  data.maxBlocks_ = 2000;
  EXPECT_NO_THROW(Eikonal::solveEikonal3D(data));
  // find the analytical solution to each vertex and compare.
  solution.resize(Eikonal::mesh_->vertices.size());
  for(size_t i = 0; i < Eikonal::mesh_->vertices.size(); i++) {
    float x = Eikonal::mesh_->vertices[i][0];
    float y = Eikonal::mesh_->vertices[i][1];
    float z = Eikonal::mesh_->vertices[i][2];
    solution[i] = std::sqrt((54. - x)*(54.-x)+(54.-y)*(54.-y)+(54.-z)*(54.-z));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  result = Eikonal::getFinalResult();
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - result[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  delete Eikonal::mesh_;
  //now test that the error decreased with each higher granularity.
  for (size_t i = 1; i < rmsError.size(); i++) {
    ASSERT_TRUE(rmsError[i - 1] >= rmsError[i]);
  }
}
