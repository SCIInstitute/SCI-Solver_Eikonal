#ifndef MESHFIM2D_EIKONAL_H
#define MESHFIM2D_EIKONAL_H

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <cusp/array1d.h>

#ifdef WIN32
#include <Windows.h>
#endif

#include <typeinfo>
#include <functional>
#include <queue>
#include <list>
#include <set>

#ifndef _EPS
#define _EPS 1e-06
#endif

typedef cusp::array1d<int, cusp::host_memory> IdxVector_h;
typedef cusp::array1d<int, cusp::device_memory> IdxVector_d;

enum LabelType2d { FarPoint = 0,ActivePoint, MaskPoint, SeedPoint,
  StopPoint, AlivePoint,ToBeAlivePoint };

class meshFIM2dEikonal {

  public:
    typedef int index;
    //  float Upwind(index vet);
    void MeshReader(char * filename);

    void SetSeedPoint(std::vector<index> SeedPoints)
    {
      m_SeedPoints = SeedPoints;
    }

    void SetMesh(TriMesh* mesh, const std::vector<float>& face_speeds, int speed_type = ONE)
    {
      m_meshPtr = mesh;
      m_meshPtr->speed_type_ = speed_type;

      orient(m_meshPtr);//  Manasi

      // have to recompute the normals and other attributes required for rendering
      if (!m_meshPtr->normals.empty()) m_meshPtr->normals.clear();//  Manasi
      m_meshPtr->need_normals();//  Manasi
      if (!m_meshPtr->adjacentfaces.empty()) m_meshPtr->adjacentfaces.clear();//  Manasi
      m_meshPtr->need_adjacentfaces();//  Manasi
      if (!m_meshPtr->across_edge.empty()) m_meshPtr->across_edge.clear();//  Manasi
      m_meshPtr->need_across_edge();//  Manasi
      if (!m_meshPtr->tstrips.empty()) m_meshPtr->tstrips.clear();//  Manasi
      m_meshPtr->need_tstrips();//  Manasi

      // initialize mesh attributes: vertT, Face.T[], Face.speedInv
      if (/*m_meshPtr->*/m_meshPtr->speed_type_ == CURVATURE)
      {
        m_meshPtr->need_curvatures();

      }
      if (/*m_meshPtr->*/m_meshPtr->speed_type_ == NOISE)
      {
        m_meshPtr->need_noise();

      }
      m_meshPtr->need_speed();

      m_meshPtr->need_faceedges();
      m_meshPtr->InitializeAttributes(face_speeds);
    }

    void InitializeLabels(int numBlock)
    {
      if (!m_meshPtr)
      {
        std::cout << "Label-vector size unknown, please set the mesh first..." << std::endl;
      }
      else
      {
        // initialize all labels to 'Far'
        size_t nv = m_meshPtr->vertices.size();
        if (m_VertLabel.size() != nv) m_VertLabel.resize(nv);
        if (m_BlockLabel.size() != numBlock) m_BlockLabel.resize(numBlock);

        for (int l = 0; l < nv; l++)
        {
          m_VertLabel[l] = FarPoint;
        }

        for (int l = 0; l < numBlock; l++)
        {
          m_BlockLabel[l] = FarPoint;
        }

        // if seeed-points are present, treat them differently
        if (!m_SeedPoints.empty())
        {
          for (int s = 0; s < m_SeedPoints.size(); s++)
          {
            m_BlockLabel[m_PartitionLabel[m_SeedPoints[s]]] = ActivePoint;//m_Label[s] = SeedPoint;
            m_VertLabel[m_SeedPoints[s]] =  SeedPoint;
            m_ActiveBlocks.insert(m_ActiveBlocks.end(), m_PartitionLabel[m_SeedPoints[s]]);
          }
        }
        else
          std::cout << "Initialize seed points before labels!!!" << std::endl;
      }
    }

    void SetStopDistance(float d)
    {
      m_StopDistance = d;

    }
    void writeVTK(std::vector< std::vector <float> > time_values);
    std::vector< std::vector< float > >  GenerateData(int numBlock,
      int maxIterations, bool verbose = false);
    void GraphPartition_METIS2(int& numBlock, int maxNumBlockVerts, bool verbose = false);
    void GraphPartition_Square(int squareLength,int squareWidth, int blockLength, int blockWidth, bool verbose = false);

    void GraphColoring();
    void PartitionFaces(int numBlock);
    meshFIM2dEikonal(){
      m_meshPtr = NULL;
    };
    ~meshFIM2dEikonal(){};

    TriMesh*                                     m_meshPtr;
    std::vector< std::set<int> >                           m_BlockNbPts;
    std::vector< std::set<int> >                           m_BlockNeighbor;
    std::vector<int>                                  m_BlockPoints;
    std::vector<int>                                  m_ColorLabel;
    int                                          m_numColor;
    std::vector<Color>                                m_faceColors;
    std::vector< std::vector<int> >                        m_PartitionFaces;
    std::vector< std::vector<int> >                        m_PartitionVerts;
    std::vector< std::vector<int> >                        m_PartitionNbFaces;
    std::vector< std::vector<TriMesh::Face> >              m_PartitionVirtualFaces;
    int                                          m_maxNumTotalFaces;
    int                                          m_maxNumVert;
    int                                          m_maxNumVertMapping;
  protected:
    std::set<index>                              m_ActiveBlocks;
    std::vector<index>                           m_SeedPoints;
    std::vector<LabelType2d>                     m_VertLabel;             // label of all the vertices active or not
    std::vector<LabelType2d>                          m_BlockLabel;            // label of blocks active or not
    float                                        m_StopDistance;
    std::vector<int> m_PartitionLabel;//TODO try to combine and use only one (npart_h on levelset)
    std::vector<int> m_BlockSizes;

    //LEVELSET vars ---- try to not have duplicates with eikonal

    int m_numPartitions;
    int m_largest_num_inside_mem;
    IdxVector_d m_xadj_d;
    IdxVector_d m_adjncy_d;
    IdxVector_d m_neighbor_sizes_d;
};
#endif
