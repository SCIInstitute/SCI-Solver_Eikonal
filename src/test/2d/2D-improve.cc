#include "gtest/gtest.h"
#include "Eikonal.h"
TEST(Improve2D, DecreaseRMSError) {
  std::vector<float> rmsError;
  std::vector< float > solution;
  float sum = 0.f;
  float radius = 100.f; //we know the radius of these spheres.
  //lowest granularity
  Eikonal* data = new Eikonal(true);
  data->filename_ = TEST_DATA_DIR + std::string("sphere_12verts.ply");
  EXPECT_NO_THROW(data->solveEikonal());
  // find the analytical solution to each vertex and compare.
  solution.resize(data->triMesh_->vertices.size());
  point first = data->triMesh_->vertices[0];
  //we know the center of these spheres.
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data->triMesh_->vertices[i][0];
    float yDot = data->triMesh_->vertices[i][1];
    float zDot = data->triMesh_->vertices[i][2];
    solution[i] = radius * std::acos(
      (first[0] * xDot + first[1] * yDot + first[2] * zDot) /
      std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot) /
      len(first));
  }
  // now calculate the RMS error for this run
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data->getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  //mid granularity
  delete data;
  data = new Eikonal(true);
  data->filename_ = TEST_DATA_DIR + std::string("sphere_74verts.ply");
  EXPECT_NO_THROW(data->solveEikonal());
  // find the analytical solution to each vertex and compare.
  solution.resize(data->triMesh_->vertices.size());
  first = data->triMesh_->vertices[0];
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data->triMesh_->vertices[i][0];
    float yDot = data->triMesh_->vertices[i][1];
    float zDot = data->triMesh_->vertices[i][2];
    solution[i] = radius * std::acos(
      (first[0] * xDot + first[1] * yDot + first[2] * zDot) /
      std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot) /
      len(first));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data->getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  //high granularity
  delete data;
  data = new Eikonal(true);
  data->filename_ = TEST_DATA_DIR + std::string("sphere_1154verts.ply");
  EXPECT_NO_THROW(data->solveEikonal());
  // find the analytical solution to each vertex and compare.
  solution.resize(data->triMesh_->vertices.size());
  first = data->triMesh_->vertices[0];
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data->triMesh_->vertices[i][0];
    float yDot = data->triMesh_->vertices[i][1];
    float zDot = data->triMesh_->vertices[i][2];
    solution[i] = radius * std::acos(
      (first[0] * xDot + first[1] * yDot + first[2] * zDot) /
      std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot) /
      len(first));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data->getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  //very high granularity
  delete data;
  data = new Eikonal(true);
  data->filename_ = TEST_DATA_DIR + std::string("sphere_4610verts.ply");
  data->maxVertsPerBlock_ = 150;
  EXPECT_NO_THROW(data->solveEikonal());
  // find the analytical solution to each vertex and compare.
  solution.resize(data->triMesh_->vertices.size());
  first = data->triMesh_->vertices[0];
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = data->triMesh_->vertices[i][0];
    float yDot = data->triMesh_->vertices[i][1];
    float zDot = data->triMesh_->vertices[i][2];
    solution[i] = radius * std::acos(
      (first[0] * xDot + first[1] * yDot + first[2] * zDot) /
      std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot) /
      len(first));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - data->getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  for (size_t i = 1; i < rmsError.size(); i++) {
    ASSERT_TRUE(rmsError[i - 1] > rmsError[i]);
  }
  delete data;
}
