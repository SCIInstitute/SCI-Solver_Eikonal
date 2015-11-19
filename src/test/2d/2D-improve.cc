#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Improve2D, DecreaseRMSError) {
  Eikonal data(true);
  std::vector<float> rmsError;
  //lowest granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_12verts.ply");
  data.maxBlocks_ = 13;
  EXPECT_NO_THROW(data.solveEikonal());
  // find the analytical solution to each vertex and compare.
  std::vector< float > solution;
  solution.resize(data.triMesh_->vertices.size());
  float radius = 19.58f; //we know the radius of these spheres.
  std::vector<float> center;
  //we know the center of these spheres.
  center.push_back(54.f);
  center.push_back(54.f);
  center.push_back(54.f);
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data.triMesh_->vertices[i][0] - center[0];
    float yDot = data.triMesh_->vertices[i][1] - center[1];
    float zDot = data.triMesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  float sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data.getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  //mid granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 11;
  EXPECT_NO_THROW(data.solveEikonal());
  // find the analytical solution to each vertex and compare.
  solution.resize(data.triMesh_->vertices.size());
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data.triMesh_->vertices[i][0] - center[0];
    float yDot = data.triMesh_->vertices[i][1] - center[1];
    float zDot = data.triMesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data.getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  //high granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_3530verts.ply");
  data.maxBlocks_ = 100;
  EXPECT_NO_THROW(data.solveEikonal());
  // find the analytical solution to each vertex and compare.
  solution.resize(data.triMesh_->vertices.size());
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data.triMesh_->vertices[i][0] - center[0];
    float yDot = data.triMesh_->vertices[i][1] - center[1];
    float zDot = data.triMesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data.getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  //very high granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_7418verts.ply");
  data.maxBlocks_ = 250;
  EXPECT_NO_THROW(data.solveEikonal());
  // find the analytical solution to each vertex and compare.
  solution.resize(data.triMesh_->vertices.size());
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data.triMesh_->vertices[i][0] - center[0];
    float yDot = data.triMesh_->vertices[i][1] - center[1];
    float zDot = data.triMesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data.getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  for (size_t i = 1; i < rmsError.size(); i++) {
    ASSERT_TRUE(rmsError[i - 1] > rmsError[i]);
  }
}
