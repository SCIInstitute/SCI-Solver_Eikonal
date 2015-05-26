#ifndef MESHFIM_H
#define MESHFIM_H


#include <tetmesh.h>
#include <typeinfo>
#include <functional>
#include <queue>
#include <list>
#include <set>
#include <time.h>
#include <stdio.h>

#ifndef _EPS
#define _EPS 1e-06
#endif


using namespace std;

enum LabelType
{
  FarPoint = 0, ActivePoint, MaskPoint, SeedPoint, StopPoint, AlivePoint, ToBeAlivePoint
};

class meshFIM
{
  public:


    typedef int index;




    //  float Upwind(index vet);
    void MeshReader(char * filename);

    //  float LocalSolver(index C, TriMesh::Face triangle);

    void SetSeedPoint(std::vector<index> SeedPoints)
    {
      m_SeedPoints = SeedPoints;

      vector<index> nb;

      if (m_meshPtr)
      {
        m_meshPtr->InitializeAttributes(0, m_SeedPoints);
        /*if (!m_SeedPoints.empty())
          {
          int ns = m_SeedPoints.size();
          for (int s = 0; s < ns; s++)
          {

          nb = m_meshPtr->neighbors[m_SeedPoints[s]];
          for (int i = 0; i<nb.size();i++)
          {
          m_ActivePoints.push_back(nb[i]);
          }


          }


          }*/

      }
    }

    void SetMesh(TetMesh* mesh)
    {
      m_meshPtr = mesh;
      m_meshPtr->vertT.resize(m_meshPtr->vertices.size());

    }

    void FindSeedPoint()
    {
      m_SeedPoints.clear();
      for (int i = 0; i < m_meshPtr->vertices.size(); i++)
      {
        double x = m_meshPtr->vertices[i][0];
        if (fabs(x - 6 ) < _EPS)
          m_SeedPoints.push_back(i);
      }
      std::cout << "SEED SIZE: " << m_SeedPoints.size() << std::endl;
    }

    void InitializePartition(int numBlock);



    void GenerateData(void);


    //void GraphPartition_Simple(int Kring, int numBlock);
    void GraphPartition_METIS2(int& numBlock, int maxNumBlockVerts);
    void GraphPartition_Square(int squareLength, int squareWidth, int squareHeight, int blockLength, int blockWidth, int blockHeight);

    void GraphColoring();
    void PartitionTets(int numBlock);

    bool gettetmem(vector<float>& tetmem, TetMesh::Tet t);
    void GetTetMem(float* &h_tetMem0, float* &h_tetMem1, float* &h_tetT);
    void GetVertMem(int* &h_vertMem, int* &h_vertMemOutside);

    void InitSpeedMat()
    {
      //    FILE* speedMatFile = fopen("speedMat6m_1.txt", "r");
      int nt = m_meshPtr->tets.size();
      //m_SpeedMat.resize(nt);
      for (int i = 0; i < nt; i++)
      {
        //        fscanf(speedMatFile, "%f %f %f %f %f %f\n", &m_meshPtr->tets[i].M[0], &m_meshPtr->tets[i].M[1], &m_meshPtr->tets[i].M[2], &m_meshPtr->tets[i].M[3], &m_meshPtr->tets[i].M[4], &m_meshPtr->tets[i].M[5]);
        m_meshPtr->tets[i].M[0] = 1;
        m_meshPtr->tets[i].M[1] = 0;
        m_meshPtr->tets[i].M[2] = 0;
        m_meshPtr->tets[i].M[3] = 1;
        m_meshPtr->tets[i].M[4] = 0;
        m_meshPtr->tets[i].M[5] = 1;

        //m_SpeedMat[i] = vector<float>(value, value + sizeof(value) / sizeof(float) );
      }

      //    fclose(speedMatFile);


    }

    void FindSeedPointLavalamp()
    {
      m_SeedPoints.clear();
      double minx = LARGENUM_TET;
      for (int i = 0; i < m_meshPtr->vertices.size(); i++)
      {
        double x = m_meshPtr->vertices[i][0];
        minx = MIN(minx, x);
        //if(abs(x+0.2) < _EPS)
        //if(fabs(x) < 0.001 && fabs(y) < 0.001 &&fabs(z) < 0.001 )
        //{
        //  m_SeedPoints.push_back(i);
        //  break;
        //}
      }
      printf("Min X is %.12f\n", minx);

      for (int i = 0; i < m_meshPtr->vertices.size(); i++)
      {
        double x = m_meshPtr->vertices[i][0];
        if (x == minx)
          m_SeedPoints.push_back(i);
      }

      if (m_SeedPoints.empty())
        printf("seed pionts list empty!!\n");
      else
        printf("Number of seed points is %lu\n", m_SeedPoints.size());

    }

    void FindSeedPointEllipse()
    {

      printf("Start FindSeedPointEllipse ...");
      m_SeedPoints.clear();
      m_SeedPointValues.clear();
      float v1, v2;
      float x1, y1, z1, x2, y2, z2;
      float K1, K2;
      K1 = 40.5;
      for (int i = 0; i < m_meshPtr->vertices.size(); i++)
      {
        //x1 = m_meshPtr->vertices[i][0] - meshsize / 2.0;
        //y1 = m_meshPtr->vertices[i][1] - meshsize / 2.0;
        //z1 = m_meshPtr->vertices[i][2] - meshsize / 2.0;
        x1 = m_meshPtr->vertices[i][0];
        y1 = m_meshPtr->vertices[i][1];
        z1 = m_meshPtr->vertices[i][2];

        v1 = x1 * x1 / (K1 * K1) + y1 * y1 * 1 / (K1 * K1) + z1 * z1 * 1 / (K1 * K1) - 1;
        //v1 = sqrt(x1*x1+ y1*y1 + z1*z1);
        //if(v1 < 3.0)
        //{
        //  m_SeedPoints.push_back(i);
        //  printf("Seed point is: vert %d.\n", i);
        //  break;
        //}

        vector< int > nbs = m_meshPtr->neighbors[i];
        for (int j = 0; j < nbs.size(); j++)
        {
          ////x2 = m_meshPtr->vertices[nbs[j]][0] - meshsize / 2.0;
          ////y2 = m_meshPtr->vertices[nbs[j]][1] - meshsize / 2.0;
          ////z2 = m_meshPtr->vertices[nbs[j]][2] - meshsize / 2.0;

          x2 = m_meshPtr->vertices[nbs[j]][0];
          y2 = m_meshPtr->vertices[nbs[j]][1];
          z2 = m_meshPtr->vertices[nbs[j]][2];


          v2 = x2 * x2 / (K1 * K1) + y2 * y2 * 1 / (K1 * K1) + z2 * z2 * 1 / (K1 * K1) - 1;

          if (v1 * v2 <= 0.0)
          {
            m_SeedPoints.push_back(i);
            K2 = sqrt(x1 * x1 + y1 * y1 * 1 + z1 * z1 * 1);
            m_SeedPointValues.push_back(fabs(K2 - K1));
            break;
          }

        }




      }

      printf("Done!\n");

    }

    meshFIM()
    {
      m_meshPtr = NULL;
    };

    ~meshFIM()
    {
    };




    TetMesh* m_meshPtr;
    vector< set<int> > m_BlockNbPts;
    vector< set<int> > m_BlockNeighbor;
    vector<int> m_BlockSizes;
    vector<int> m_BlockPoints;
    vector<int> m_ColorLabel;
    int m_numColor;
    vector< vector<int> > m_PartitionTets;
    //vector< vector<int> >                        m_PartitionVerts;
    vector< vector<int> > m_PartitionInVerts;
    vector< set<int> > m_PartitionOutVerts;
    vector< vector<int> > m_PartitionNbTets;
    vector< vector< vector<int> > > m_VertMapping;
    vector< vector<int> > m_OutVertMapping;
    vector< vector<int> > m_OneRingStrip;
    vector< vector< vector<float> > > m_tetMem;
    vector< vector<int> > m_blockVertMapping;
    vector< vector<int> > m_blockVertMappingInside;
    vector< vector<int> > m_blockVertMappingOutside;

    vector< vector<TetMesh::Tet> > m_PartitionVirtualTets;
    int m_maxNumTotalTets;
    //int                                          m_maxNumVert;
    int m_maxNumInVert;
    int m_maxNumVertMapping;
    int m_maxVertMappingInside;
    int m_maxVertMappingOutside;
    int m_maxNumOutVertMapping;
    int m_maxLengthStrip;
    int m_numBlock;
    vector<float> m_SeedPointValues;

    std::vector<int> m_PartitionLabel; // label of vertices belong to which partition

  protected:

    std::set<index> m_ActiveBlocks;
    std::vector<index> m_SeedPoints;
    std::vector<LabelType> m_VertLabel; // label of all the vertices active or not
    vector<LabelType> m_BlockLabel; // label of blocks active or not
    float m_StopDistance;
};

#endif
