/*
   Szymon Rusinkiewicz
   Princeton University

   TriMesh_connectivity.cc
   Manipulate data structures that describe connectivity between faces and verts.
 */


#include <stdio.h>
#include <time.h>
#include "TriMesh.h"
#include <algorithm>

void TriMesh::need_Rinscribe()
{
  need_faceedges();
  int nf = faces.size();
  if (!radiusInscribe.empty())
  {
    return;
  }

  radiusInscribe.resize(nf);

  for (int i=0; i<nf; i++)
  {
    Face f = faces[i];
    double e1 = f.edgeLens[0];
    double e2 = f.edgeLens[1];
    double e3 = f.edgeLens[2];
    double s = (e1+e2+e3)/2.0;
    radiusInscribe[i] = sqrt(s *(s-e1) * (s-e2) * (s - e3)) / s;
  }

}


void TriMesh::need_faces()
{
  if (!faces.empty())
  {
    return;
  }
  if (!tstrips.empty())
    unpack_tstrips();
  else if (!grid.empty())
    triangulate_grid();


}


void TriMesh::need_faceedges()
{
  if (faces.empty())
  {
    printf("No faces to compute face edges!!!\n");
    return;
  }
  int numFaces = faces.size();
  for (int i = 0; i < numFaces; i++)
  {
    Face f = faces[i];
    point edge01 = vertices[f[1]] - vertices[f[0]];
    point edge12 = vertices[f[2]] - vertices[f[1]];
    point edge20 = vertices[f[0]] - vertices[f[2]];
    faces[i].edgeLens[0] =sqrt(edge01[0]*edge01[0] + edge01[1]*edge01[1] + edge01[2]*edge01[2]);
    faces[i].edgeLens[1] =sqrt(edge12[0]*edge12[0] + edge12[1]*edge12[1] + edge12[2]*edge12[2]);
    faces[i].edgeLens[2] =sqrt(edge20[0]*edge20[0] + edge20[1]*edge20[1] + edge20[2]*edge20[2]);


  }

}

void TriMesh::need_noise(float low, float high)
{
  noiseOnVert.clear();
  need_neighbors();
  int nv = vertices.size();
  noiseOnVert.resize(nv);
  srand( (unsigned)time( NULL ) );

  for (int i = 0;i<nv; i++) {
    noiseOnVert[i] = (float)rand() /
      (RAND_MAX)*(high - low) + low;  //random number between [low,high]
  }

  //iterate
  int iterNum = 0;
  for (int i=0;i<iterNum; i++) {
    for (int j=0; j<nv;j++) {
      noiseOnVert[j] = 0;
      std::vector<int> nb = neighbors[j];
      for (int k=0;k<nb.size();k++) {
        noiseOnVert[j] +=noiseOnVert[neighbors[j][k]];
      }
      noiseOnVert[j] /= nb.size();
    }
  }
}

void TriMesh::need_speed()
{
  //this is a set of optional speed functions to use.
  int nf = faces.size();

  for (int i =0; i<nf;i++)
  {
    switch (this->speed_type_)
    {
    case CURVATURE:
      faces[i].speedInv = (abs(curv1[faces[i][0]] +
            curv2[faces[i][0]]) +
          std::abs(curv1[faces[i][1]] + curv2[faces[i][1]]) +
          std::abs(curv1[faces[i][2]] +
            curv2[faces[i][2]])) / 6.0;
      break;
    case ONE:
      faces[i].speedInv = 1.0;
      break;
    case NOISE:
      faces[i].speedInv =( noiseOnVert[faces[i][0]] +
          noiseOnVert[faces[i][1]] +
          noiseOnVert[faces[i][2]] ) / 3;
      break;
    default:
      faces[i].speedInv = 1.0;
      break;
    }
  }
}


// Find the direct neighbors of each vertex
void TriMesh::need_neighbors(bool verbose)
{
  if (!neighbors.empty())
    return;
  need_faces();

  if (verbose)
    printf("Finding vertex neighbors... ");
  int nv = vertices.size(), nf = faces.size();

  std::vector<int> numneighbors(nv);
  for (int i = 0; i < nf; i++) {
    numneighbors[faces[i][0]]++;
    numneighbors[faces[i][1]]++;
    numneighbors[faces[i][2]]++;
  }

  neighbors.resize(nv);
  for (int i = 0; i < nv; i++)
    neighbors[i].reserve(numneighbors[i]+2); // Slop for boundaries

  for (int i = 0; i < nf; i++) {
    for (int j = 0; j < 3; j++) {
      std::vector<int> &me = neighbors[faces[i][j]];
      int n1 = faces[i][(j+1)%3];
      int n2 = faces[i][(j+2)%3];
      if (find(me.begin(), me.end(), n1) == me.end())
        me.push_back(n1);
      if (find(me.begin(), me.end(), n2) == me.end())
        me.push_back(n2);
    }
  }

  if(verbose)
    printf("Done.\n");
}


// Find the faces touching each vertex
void TriMesh::need_adjacentfaces(bool verbose)
{
  if (!adjacentfaces.empty())
    return;
  need_faces();

  if (verbose)
    printf("Finding vertex to triangle maps... ");
  int nv = vertices.size(), nf = faces.size();

  std::vector<int> numadjacentfaces(nv);
  for (int i = 0; i < nf; i++) {
    numadjacentfaces[faces[i][0]]++;
    numadjacentfaces[faces[i][1]]++;
    numadjacentfaces[faces[i][2]]++;
  }

  adjacentfaces.resize(vertices.size());
  for (int i = 0; i < nv; i++)
    adjacentfaces[i].reserve(numadjacentfaces[i]);

  for (int i = 0; i < nf; i++) {
    for (int j = 0; j < 3; j++)
      adjacentfaces[faces[i][j]].push_back(i);
  }

  if(verbose)
    printf("Done.\n");
}

void TriMesh::need_face_virtual_faces()
{

  std::vector<Face> t_faces;
  Face f;
  int numFaces = faces.size();
  faceVirtualFaces.resize(numFaces);
  for (int i = 0; i < numFaces; i++)
  {
    t_faces.clear();
    f = faces[i];

    for (int j = 0; j< 3 ; j++)
    {
      if(!IsNonObtuse(f[j],f))// check angle: if non-obtuse, return existing face
      {
        int nfae = across_edge[i][j];
        if (nfae > -1)
        {
          SplitFace(t_faces,f[j],f,nfae);// if obtuse, split face till we get all acute angles
        }
        else
          printf("NO cross edge!!! Maybe a hole!!\n");

      }
    }

    faceVirtualFaces[i] = t_faces;
  }
}


// Find the face across each edge from each other face (-1 on boundary)
// If topology is bad, not necessarily what one would expect...
void TriMesh::need_across_edge()
{
  if (!across_edge.empty())
    return;
  need_adjacentfaces();

  dprintf("Finding across-edge maps... ");

  int nf = faces.size();
  across_edge.resize(nf, Face(-1,-1,-1));

  for (int i = 0; i < nf; i++) {
    for (int j = 0; j < 3; j++) {
      if (across_edge[i][j] != -1)
        continue;
      int v1 = faces[i][(j+1)%3];
      int v2 = faces[i][(j+2)%3];
      const std::vector<int> &a1 = adjacentfaces[v1];
      const std::vector<int> &a2 = adjacentfaces[v2];
      for (int k1 = 0; k1 < a1.size(); k1++) {
        int other = a1[k1];
        if (other == i)
          continue;
        std::vector<int>::const_iterator it =
          std::find(a2.begin(), a2.end(), other);
        if (it == a2.end())
          continue;
        int ind = (faces[other].indexof(v1)+1)%3;
        if (faces[other][(ind+1)%3] != v2)
          continue;
        across_edge[i][j] = other;
        across_edge[other][ind] = i;
        break;
      }
    }
  }

  dprintf("Done.\n");
}

