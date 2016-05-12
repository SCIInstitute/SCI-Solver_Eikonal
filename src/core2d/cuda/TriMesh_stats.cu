/*
   Szymon Rusinkiewicz
   Princeton University

   TriMesh_stats.cc
   Computation of various statistics on the mesh.
 */

#include "TriMesh.h"
#include <algorithm>

// A characteristic "feature size" for the mesh.  Computed as an approximation
// to the median edge length
float TriMesh::feature_size()
{
  need_faces();
  if (faces.empty())
    return 0.0f;

  int nf = faces.size();
  int nsamp = std::min(nf / 2, 333);

  std::vector<float> samples;
  samples.reserve(nsamp * 3);

  for (int i = 0; i < nsamp; i++) {
    // Quick 'n dirty portable random number generator
    static unsigned randq = 0;
    randq = unsigned(1664525) * randq + unsigned(1013904223);

    int ind = randq % nf;
    const point &p0 = vertices[faces[ind][0]];
    const point &p1 = vertices[faces[ind][1]];
    const point &p2 = vertices[faces[ind][2]];
    samples.push_back(dist2(p0,p1));
    samples.push_back(dist2(p1,p2));
    samples.push_back(dist2(p2,p0));
  }
  nth_element(samples.begin(),
      samples.begin() + samples.size()/2,
      samples.end());
  return sqrt(samples[samples.size()/2]);
}


void TriMesh::InitializeAttributes(const std::vector<float>& face_speeds)
{
  // pre-compute faces, normals, and other per-vertex properties that may be needed
  this->need_neighbors();
  this->need_normals();
  this->need_adjacentfaces();
  this->need_across_edge();
  this->need_faces();
  this->need_face_virtual_faces();

  // for all faces: initialize per-vertex travel time and face-speed
  size_t nf = this->faces.size();
  for (int f = 0; f < nf; f++)
  {
    // travel time
    this->faces[f].T[0] = this->vertT[this->faces[f][0]];
    this->faces[f].T[1] = this->vertT[this->faces[f][1]];
    this->faces[f].T[2] = this->vertT[this->faces[f][2]];
    if (face_speeds.empty() || face_speeds.size() != this->faces.size()){
      this->faces[f].speedInv = 1.f;
    } else {
      this->faces[f].speedInv = 1.f / face_speeds[f];
    }
  }
}
