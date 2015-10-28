#include "gtest/gtest.h"
#include "Eikonal2D.h"
TEST(Improve2D, DecreaseRMSError) {
  Eikonal2D::Eikonal2D data;
  std::vector<float> rmsError;
  //lowest granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_12verts.ply");
  data.maxBlocks_ = 2;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  // find the analytical solution to each vertex and compare.
  std::vector< float > solution;
  solution.resize(Eikonal2D::mesh_->vertices.size());
  float radius = 19.58f; //we know the radius of these spheres.
  std::vector<float> center;
  //we know the center of these spheres.
  center.push_back(54.f);
  center.push_back(54.f);
  center.push_back(54.f);
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = Eikonal2D::mesh_->vertices[i][0] - center[0];
    float yDot = Eikonal2D::mesh_->vertices[i][1] - center[1];
    float zDot = Eikonal2D::mesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  float sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - Eikonal2D::getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  delete Eikonal2D::mesh_;
  //mid granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_266verts.ply");
  data.maxBlocks_ = 11;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  // find the analytical solution to each vertex and compare.
  solution.resize(Eikonal2D::mesh_->vertices.size());
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = Eikonal2D::mesh_->vertices[i][0] - center[0];
    float yDot = Eikonal2D::mesh_->vertices[i][1] - center[1];
    float zDot = Eikonal2D::mesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - Eikonal2D::getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  delete Eikonal2D::mesh_;
  //high granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_3530verts.ply");
  data.maxBlocks_ = 100;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  // find the analytical solution to each vertex and compare.
  solution.resize(Eikonal2D::mesh_->vertices.size());
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = Eikonal2D::mesh_->vertices[i][0] - center[0];
    float yDot = Eikonal2D::mesh_->vertices[i][1] - center[1];
    float zDot = Eikonal2D::mesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - Eikonal2D::getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  delete Eikonal2D::mesh_;
  //very high granularity
  data.filename_ = TEST_DATA_DIR + std::string("sphere_7418verts.ply");
  data.maxBlocks_ = 250;
  EXPECT_NO_THROW(Eikonal2D::solveEikonal2D(data));
  // find the analytical solution to each vertex and compare.
  solution.resize(Eikonal2D::mesh_->vertices.size());
  for (size_t i = 0; i < solution.size(); i++) {
    float xDot = Eikonal2D::mesh_->vertices[i][0] - center[0];
    float yDot = Eikonal2D::mesh_->vertices[i][1] - center[1];
    float zDot = Eikonal2D::mesh_->vertices[i][2] - center[2];
    solution[i] = radius * std::acos( zDot /
        std::sqrt(xDot * xDot + yDot * yDot + zDot * zDot));
  }
  // now calculate the RMS error for this run
  sum = 0.f;
  for (size_t j = 0; j < solution.size(); j++) {
    float err = solution[j] - Eikonal2D::getFinalResult()[j];
    sum +=  err * err;
  }
  rmsError.push_back(std::sqrt(sum / static_cast<float>(solution.size())));
  delete Eikonal2D::mesh_;
  for (size_t i = 1; i < rmsError.size(); i++) {
    ASSERT_TRUE(rmsError[i - 1] > rmsError[i]);
  }
}
