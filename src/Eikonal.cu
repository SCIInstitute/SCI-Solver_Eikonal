#include "Eikonal.h"
#include <cmath>

Eikonal::Eikonal(bool isTriMesh, std::string fname, bool verbose) :
  verbose_(verbose),
  filename_(fname),
  maxBlocks_(16),
  maxVertsPerBlock_(64),
  stopDistance_(50000.f),
  isStructured_(false),
  userSetInitial_(false),
  speedType_(ONE),
  squareLength_(16),
  squareWidth_(16),
  squareDepth_(16),
  squareBlockLength_(1),
  squareBlockWidth_(1),
  squareBlockDepth_(1),
  maxIterations_(100),
  triMesh_(NULL),
  tetMesh_(NULL),
  FIMPtr2d_(NULL),
  FIMPtr3d_(NULL),
  isTriMesh_(isTriMesh) {}

  Eikonal::~Eikonal() {
    if (this->tetMesh_ != NULL)
      delete this->tetMesh_;
    if (this->triMesh_ != NULL)
      delete this->triMesh_;
    if (this->FIMPtr2d_ != NULL)
      delete this->FIMPtr2d_;
    if (this->FIMPtr3d_ != NULL)
      delete this->FIMPtr3d_;
  }


std::vector < float >&  Eikonal::getFinalResult() {
  return iteration_values_.at(iteration_values_.size() - 1);
}

std::vector < float >&  Eikonal::getResultAtIteration(size_t i) {
  return iteration_values_.at(i);
}

size_t  Eikonal::numIterations() {
  return iteration_values_.size();
}

void Eikonal::writeVTK(bool all) {
  std::vector<std::vector<float> > vals;
  for(size_t i = all?0:this->iteration_values_.size()  - 1;
      i < this->iteration_values_.size(); i++) {
    vals.push_back(this->iteration_values_[i]);
  }
  if (FIMPtr2d_ != NULL)
    FIMPtr2d_->writeVTK(vals);
  else
    FIMPtr3d_->writeVTK(vals);
}

void Eikonal::initSpeedMtxMultipliers(
    std::vector<float> values) {
  // num triangles size for 2D
  // num tets or num tets * 6 size for 3D
  this->speedMtxMultipliers_ = values;
}

void Eikonal::initializeVertices(std::vector<float> values) {
  if (this->triMesh_ == NULL && this->tetMesh_ == NULL) {
    std::cerr << "You must initialize the mesh first!" << std::endl;
    exit(0);
  }
  if (this->triMesh_ != NULL) {
    if (values.size() != this->triMesh_->vertices.size()) {
      std::cerr << "Initialize values size does not match number of vertices!"
        << std::endl;
      exit(0);
    }
    this->triMesh_->vertT.resize(this->triMesh_->vertices.size());
    for (size_t i = 0; i < values.size(); i++) {
      this->triMesh_->vertT[i] = values[i];
    }
  } else {
    if (values.size() != this->tetMesh_->vertices.size()) {
      std::cerr << "Initialize values size does not match number of vertices!"
        << std::endl;
      exit(0);
    }
    this->tetMesh_->vertT.resize(this->tetMesh_->vertices.size());
    for (size_t i = 0; i < values.size(); i++) {
      this->tetMesh_->vertT[i] = values[i];
    }

  }
  this->userSetInitial_ = true;
}

void Eikonal::initializeMesh() {
  if (this->isTriMesh_) {
    if (this->triMesh_ == NULL) {
      this->triMesh_ = TriMesh::read(this->filename_.c_str(), this->verbose_);
      if (this->triMesh_ == NULL)
      {
        printf("File open failed!!\n");
        exit(1);
      }
    }
  } else {
    if (this->tetMesh_ == NULL) {
      tetgenio in;
      if (!(in.load_tetmesh((char*)this->filename_.c_str(),
              this->verbose_))) {
        exit(1);
      }

      this->tetMesh_ = new TetMesh();
      this->tetMesh_->init(
          in.pointlist,
          in.numberofpoints,
          in.trifacelist,
          in.numberoffacets,
          in.tetrahedronlist,
          in.numberoftetrahedra,
          in.tetrahedronattributelist,
          this->speedMtxMultipliers_,
          this->verbose_);
      this->tetMesh_->need_neighbors(this->verbose_);
      this->tetMesh_->need_adjacenttets(this->verbose_);
      this->tetMesh_->need_tet_virtual_tets(this->verbose_);
    }
  }
}

void Eikonal::solveEikonal() {
  clock_t starttime, endtime;
  starttime = clock();
  if (this->isTriMesh_) {
    if (this->triMesh_ == NULL) {
      this->initializeMesh();
    }
    FIMPtr2d_ = new meshFIM2dEikonal;
    //initialize the first point as the "Seed"
    if (!this->userSetInitial_) {
      this->triMesh_->vertT.resize(this->triMesh_->vertices.size());
      this->triMesh_->vertT[0] = 0.;
      for (size_t i = 1; i < this->triMesh_->vertices.size(); i++) {
        this->triMesh_->vertT[i] = LARGENUM;
      }
      FIMPtr2d_->SetSeedPoint(std::vector<int>(1, 0));
    } else {
      std::vector<int> found_seeds;
      for (size_t i = 0; i < this->triMesh_->vertices.size(); i++) {
        if (this->triMesh_->vertT[i] == 0.) {
          found_seeds.push_back(static_cast<int>(i));
        }
      }
      FIMPtr2d_->SetSeedPoint(found_seeds);
    }
    FIMPtr2d_->SetMesh(this->triMesh_, this->speedMtxMultipliers_, this->speedType_);
    FIMPtr2d_->SetStopDistance(this->stopDistance_);
    if (this->isStructured_) {
      int numBlockLength = (this->squareLength_ / this->squareBlockLength_);
      int numBlockWidth = (this->squareWidth_ / this->squareBlockWidth_);
      this->maxBlocks_ = numBlockLength * numBlockWidth;
      FIMPtr2d_->GraphPartition_Square(this->squareLength_, this->squareWidth_,
          this->squareBlockLength_,
          this->squareBlockWidth_,
          this->verbose_);

    } else {
      FIMPtr2d_->GraphPartition_METIS2(this->maxBlocks_,
          this->maxVertsPerBlock_, this->verbose_);
    }
    FIMPtr2d_->PartitionFaces(this->maxBlocks_);
    FIMPtr2d_->InitializeLabels(this->maxBlocks_);

    iteration_values_ =
      FIMPtr2d_->GenerateData(this->maxBlocks_, this->maxIterations_, this->verbose_);
  } else {
    if (this->tetMesh_ == NULL) {
      this->initializeMesh();
    }
    FIMPtr3d_ = new meshFIM3dEikonal;
    FIMPtr3d_->SetMesh(this->tetMesh_);
    //initialize the first point as the "Seed"
    if (!this->userSetInitial_) {
      this->tetMesh_->vertT.resize(this->tetMesh_->vertices.size());
      this->tetMesh_->vertT[0] = 0.;
      for (size_t i = 1; i < this->tetMesh_->vertices.size(); i++) {
        this->tetMesh_->vertT[i] = LARGENUM;
      }
      FIMPtr3d_->SetSeedPoint(std::vector<int>(1, 0));
    } else {
      std::vector<int> found_seeds;
      for (size_t i = 0; i < this->tetMesh_->vertices.size(); i++) {
        if (this->tetMesh_->vertT[i] == 0.) {
          found_seeds.push_back(static_cast<int>(i));
        }
      }
      FIMPtr3d_->SetSeedPoint(found_seeds);
    }

    if (this->isStructured_) {
      int numBlockLength = (this->squareLength_ / this->squareBlockLength_);
      int numBlockWidth = (this->squareWidth_ / this->squareBlockWidth_);
      int numBlockDepth  = (this->squareDepth_ / this->squareBlockDepth_);
      this->maxBlocks_ = numBlockLength * numBlockWidth * numBlockDepth;
      FIMPtr3d_->GraphPartition_Square(this->squareLength_,
          this->squareWidth_,
          this->squareDepth_,
          this->squareBlockLength_,
          this->squareBlockWidth_,
          this->squareBlockDepth_,
          this->verbose_);

    } else {
      FIMPtr3d_->GraphPartition_METIS2(this->maxBlocks_,
          this->maxVertsPerBlock_, this->verbose_);
    }
    FIMPtr3d_->m_numBlock = this->maxBlocks_;
    FIMPtr3d_->PartitionTets(this->maxBlocks_, this->verbose_);
    iteration_values_ =
      FIMPtr3d_->GenerateData(this->maxIterations_, this->verbose_);
  }
  endtime = clock();
  double duration = (double)(endtime - starttime) * 1000 / CLOCKS_PER_SEC;

  if (this->verbose_)
    printf("Computing time : %.10lf ms\n", duration);
}

void Eikonal::printErrorGraph(std::vector<float> solution) {
  // now calculate the RMS error for each iteration
  std::vector<float> rmsError;
  rmsError.resize(numIterations());
  for (size_t i = 0; i < numIterations(); i++) {
    float sum = 0.f;
    std::vector<float> result = getResultAtIteration(i);
    for (size_t j = 0; j < solution.size(); j++) {
      float err = std::abs(solution[j] - result[j]);
      sum += err * err;
    }
    rmsError[i] = std::sqrt(sum / static_cast<float>(solution.size()));
  }
  //determine the log range
  float max_err = rmsError[0];
  float min_err = rmsError[rmsError.size() - 1];
  int max_log = -10, min_log = 10;
  while (std::pow(static_cast<float>(10),
        static_cast<float>(max_log)) < max_err) max_log++;
  while (std::pow(static_cast<float>(10),
        static_cast<float>(min_log)) > min_err) min_log--;
  // print the error graph
  printf("\n\nlog(Err)|\n");
  bool printTick = true;
  for (int i = max_log; i >= min_log; i--) {
    if (printTick) {
      printf("   10^%2d|", i);
    } else {
      printf("        |");
    }
    for (size_t j = 0; j < numIterations(); j++) {
      if (rmsError[j] > std::pow(static_cast<float>(10),
            static_cast<float>(i)) &&
          rmsError[j] < std::pow(static_cast<float>(10),
            static_cast<float>(i + 1)))
        printf("*");
      else
        printf(" ");
    }
    printf("\n");
    printTick = !printTick;
  }
  printf("--------|------------------------------------------");
  printf("  Converged to: %.4f\n", rmsError[rmsError.size() - 1]);
  printf("        |1   5    10   15   20   25   30   35\n");
  printf("                   Iteration\n");
}
