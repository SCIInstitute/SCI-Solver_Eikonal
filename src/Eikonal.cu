#include "Eikonal.h"

Eikonal::Eikonal(bool isTriMesh, std::string fname, bool verbose) :
  verbose_(verbose),
  filename_(fname),
  seedPointList_(std::vector<int>(1, 0)),
  maxBlocks_(10003),
  maxVertsPerBlock_(64),
  stopDistance_(50000.f),
  isStructured_(false),
  speedType_(ONE),
  squareLength_(16),
  squareWidth_(16),
  squareDepth_(16),
  squareBlockLength_(1),
  squareBlockWidth_(1),
  squareBlockDepth_(1),
  maxIterations_(1000),
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

void Eikonal::writeVTK() {
  if (FIMPtr2d_ != NULL)
    FIMPtr2d_->writeVTK(this->iteration_values_);
  else
    FIMPtr3d_->writeVTK(this->iteration_values_);
}


void Eikonal::solveEikonal() {
  clock_t starttime, endtime;
  starttime = clock();
  if (this->isTriMesh_) {
    this->triMesh_ = TriMesh::read(this->filename_.c_str(), this->verbose_);
    if (!this->triMesh_) exit(1);
    FIMPtr2d_ = new meshFIM2dEikonal;
    FIMPtr2d_->SetSeedPoint(this->seedPointList_);
    FIMPtr2d_->SetMesh(this->triMesh_, this->speedType_);
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
        in.numberoftetrahedronattributes,
        in.tetrahedronattributelist, this->verbose_);
    this->tetMesh_->need_neighbors(this->verbose_);
    this->tetMesh_->need_adjacenttets(this->verbose_);
    this->tetMesh_->need_tet_virtual_tets(this->verbose_);

    FIMPtr3d_ = new meshFIM3dEikonal;
    FIMPtr3d_->SetSeedPoint(this->seedPointList_);
    FIMPtr3d_->SetMesh(this->tetMesh_);
    FIMPtr3d_->InitSpeedMat();

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
  while (std::pow(static_cast<float>(10), max_log) < max_err) max_log++;
  while (std::pow(static_cast<float>(10), min_log) > min_err) min_log--;
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
      if (rmsError[j] > std::pow(static_cast<float>(10), i) &&
          rmsError[j] < std::pow(static_cast<float>(10), i + 1))
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