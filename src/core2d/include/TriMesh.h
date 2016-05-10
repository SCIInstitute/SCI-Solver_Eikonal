#ifndef TRIMESH_H
#define TRIMESH_H
/*
   Szymon Rusinkiewicz
   Princeton University

   TriMesh.h
   Class for triangle meshes.
 */

#ifndef LARGENUM
#define LARGENUM 10000000.0
#endif

#ifndef ONE
#define ONE 1
#endif
#ifndef CURVATURE
#define CURVATURE 2
#endif
#ifndef NOISE
#define NOISE 3
#endif
#ifndef MANUAL_SPEED
#define MANUAL_SPEED 3
#endif

#include <cstdio>
#include "Vec.h"
#include "Color.h"
#include "math.h"
#include <vector>
#include <list>

#ifndef M_PI
#define M_PI 3.14159265359
#endif
class TriMesh {
  protected:
    static bool read_helper(const char *filename, TriMesh *mesh);

  public:
    // Types
    struct Face {
      int v[3];
      float speedInv;
      float T[3];
      Vec<3,float> edgeLens;  // edge length for 01, 12, 20
      int material_;
      Face() : speedInv(1.f), material_(0) {}
      Face(const int &v0, const int &v1, const int &v2) : speedInv(1.f), material_(0)
      {
        v[0] = v0; v[1] = v1; v[2] = v2;
      }
      Face(const int *v_) : speedInv(1.f), material_(0)
      {
        v[0] = v_[0]; v[1] = v_[1]; v[2] = v_[2];
      }
      int &operator[] (int i) { return v[i]; }
      const int &operator[] (int i) const { return v[i]; }
      operator const int * () const { return &(v[0]); }
      operator const int * () { return &(v[0]); }
      operator int * () { return &(v[0]); }
      int indexof(int v_) const
      {
        return (v[0] == v_) ? 0 :
          (v[1] == v_) ? 1 :
          (v[2] == v_) ? 2 : -1;
      }
    };

    struct BBox {
      point min, max;
      point center() const { return (min+max) * 0.5f; }
      vec size() const { return max - min; }
      bool valid;
      BBox() : valid(false)
      {}
    };

    struct BSphere {
      point center;
      float r;
      bool valid;
      BSphere() : valid(false)
      {}
    };

    // Enums
    enum tstrip_rep { TSTRIP_LENGTH, TSTRIP_TERM };
    enum { GRID_INVALID = -1 };

    // The basics: vertices and faces
    std::vector<point> vertices;
    std::vector<Face> faces;

    // Triangle strips
    std::vector<int> tstrips;

    // Grid, if present
    std::vector<int> grid;
    int grid_width, grid_height;

    // Other per-vertex properties
    std::vector<Color> colors;
    std::vector<float> confidences;
    std::vector<unsigned> flags;
    unsigned flag_curr;

    // Computed per-vertex properties
    std::vector<vec> normals;
    std::vector<vec> pdir1, pdir2;
    std::vector<float> curv1, curv2;
    std::vector< Vec<4,float> > dcurv;
    std::vector<vec> cornerareas;
    std::vector<float> pointareas;
    std::vector<float> vertT;

    //the type of speed on faces
    int speed_type_;

    // Bounding structures
    BBox bbox;
    BSphere bsphere;

    std::vector<float> noiseOnVert;

    std::vector< std::vector<Face> > faceVirtualFaces;

    // Connectivity structures:
    //  For each vertex, all neighboring vertices
    std::vector< std::vector<int> > neighbors;

    std::vector< std::vector<int*> > NonObtuseNeighborFaces;

    //  For each vertex, all neighboring faces
    std::vector< std::vector<int> > adjacentfaces;
    std::vector<float> radiusInscribe;
    //  For each face, the three faces attached to its edges
    //  (for example, across_edge[3][2] is the number of the face
    //   that's touching the edge opposite vertex 2 of face 3)
    std::vector<Face> across_edge;

    // Compute all this stuff...
    void need_tstrips();
    void convert_strips(tstrip_rep rep);
    void unpack_tstrips();
    void need_noise(float low = 0.f, float high = 1.f);
    void need_speed();
    void triangulate_grid();
    void need_faces();
    void need_faceedges();
    void need_normals();
    void need_pointareas();
    void need_curvatures();
    void need_dcurv();
    void need_bbox();
    void need_bsphere();
    void need_neighbors(bool verbose = false);
    void need_adjacentfaces(bool verbose = false);
    void need_Rinscribe();
    void need_across_edge();
    void need_face_virtual_faces();

    // Input and output
    static TriMesh *read(const char *filename, bool verbose = false);
    void write(const char *filename);

    // Statistics
    // XXX - Add stuff here
    float feature_size();

    // Useful queries
    // XXX - Add stuff here
    bool is_bdy(int v)
    {
      if (neighbors.empty()) need_neighbors();
      if (adjacentfaces.empty()) need_adjacentfaces();
      return neighbors[v].size() != adjacentfaces[v].size();
    }


    vec trinorm(int f)
    {
      if (faces.empty()) need_faces();
      return ::trinorm(vertices[faces[f][0]], vertices[faces[f][1]],
          vertices[faces[f][2]]);
    }

    // FIM: check angle for at a given vertex, for a given face
    bool IsNonObtuse(int v, Face f)
    {
      int iV = f.indexof(v);

      point A = this->vertices[v];
      point B = this->vertices[f[(iV+1)%3]];
      point C = this->vertices[f[(iV+2)%3]];

      float a = dist(B,C);
      float b = dist(A,C);
      float c = dist(A,B);

      float angA = 0.0;//acos( (b*b + c*c - a*a) / (2*b*c) );   //manasi

      if ((a > 0) && (b > 0) && (c > 0))//  Manasi stack overflow
      {//  Manasi stack overflow
        angA = acos( (b*b + c*c - a*a) / (2*b*c) );//  Manasi stack overflow
      }//  Manasi stack overflow

      return (angA < M_PI/2.0f);
      //return 1;
    }

    // FIM: given a vertex, find an all-acute neighborhood of faces
    void SplitFace(std::vector<Face> &acFaces, int v, Face cf, int nfAdj)
    {
      // get all the four vertices in order
      /* v1         v4
         +-------+
         \     . \
         \   .   \
         \ .     \
         +-------+
         v2         v3 */

      int iV = cf.indexof(v);                      // get index of v in terms of cf
      int v1 = v;
      int v2 = cf[(iV+1)%3];
      int v4 = cf[(iV+2)%3];
      iV = this->faces[nfAdj].indexof(v2);        // get index of v in terms of adjacent face
      int v3 = this->faces[nfAdj][(iV+1)%3];

      // create faces (v1,v3,v4) and (v1,v2,v3), check angle at v1
      Face f1(v1,v3,v4);

      Face f2(v1,v2,v3);

      if (IsNonObtuse(v,f1))
      {
        switch (this->speed_type_)
        {

        case CURVATURE:
          f1.speedInv = ( std::abs(curv1[f1[0]] + curv2[f1[0]]) + 
            std::abs(curv1[f1[1]] + curv2[f1[1]]) + 
            std::abs(curv1[f1[2]] + curv2[f1[2]]) ) / 6;
          break;
        case ONE:
          f1.speedInv = 1.0;
          break;
        case NOISE:
          f1.speedInv =( noiseOnVert[f1[0]] + noiseOnVert[f1[1]] + noiseOnVert[f1[2]] )/ 3;
          break;
        default:
          f1.speedInv = cf.speedInv;
          break;

        }

        point edge01 = vertices[f1[1]] - vertices[f1[0]];
        point edge12 = vertices[f1[2]] - vertices[f1[1]];
        point edge20 = vertices[f1[0]] - vertices[f1[2]];
        f1.edgeLens[0] =sqrt(edge01[0]*edge01[0] + edge01[1]*edge01[1] + edge01[2]*edge01[2]);
        f1.edgeLens[1] =sqrt(edge12[0]*edge12[0] + edge12[1]*edge12[1] + edge12[2]*edge12[2]);
        f1.edgeLens[2] =sqrt(edge20[0]*edge20[0] + edge20[1]*edge20[1] + edge20[2]*edge20[2]);
        acFaces.push_back(f1);
      }
      else
      {
        int nfAdj_new = this->across_edge[nfAdj][this->faces[nfAdj].indexof(v2)];
        if (nfAdj_new > -1)
          SplitFace(acFaces,v,f1,nfAdj_new);
        else
          printf("NO cross edge!!! Maybe a hole!!\n");


      }

      if (IsNonObtuse(v,f2))
      {
        switch (this->speed_type_)
        {
        case CURVATURE:
          f2.speedInv = ( std::abs(curv1[f2[0]] + curv2[f2[0]]) + std::abs( curv1[f2[1]] + curv2[f2[1]]) + std::abs(curv1[f2[2]] + curv2[f2[2]]) ) / 6;
          break;
        case ONE:
          f2.speedInv = 1.0;
          break;
        case NOISE:
          f2.speedInv = ( noiseOnVert[f2[0]] + noiseOnVert[f2[1]] + noiseOnVert[f2[2]] )/ 3;
          break;
        default:
          f2.speedInv = cf.speedInv;
          break;

        }
        point edge01 = vertices[f2[1]] - vertices[f2[0]];
        point edge12 = vertices[f2[2]] - vertices[f2[1]];
        point edge20 = vertices[f2[0]] - vertices[f2[2]];
        f2.edgeLens[0] =sqrt(edge01[0]*edge01[0] + edge01[1]*edge01[1] + edge01[2]*edge01[2]);
        f2.edgeLens[1] =sqrt(edge12[0]*edge12[0] + edge12[1]*edge12[1] + edge12[2]*edge12[2]);
        f2.edgeLens[2] =sqrt(edge20[0]*edge20[0] + edge20[1]*edge20[1] + edge20[2]*edge20[2]);
        acFaces.push_back(f2);
      }
      else
      {
        int nfAdj_new = this->across_edge[nfAdj][this->faces[nfAdj].indexof(v4)];
        if(nfAdj_new > -1)
          SplitFace(acFaces,v,f2,nfAdj_new);
        else
          printf("NO cross edge!!! Maybe a hole!!\n");
      }
    }

    // FIM: one ring function
    std::vector<Face> GetOneRing(int v)
    {
      // make sure we have the across-edge map
      if (this->across_edge.empty())
        this->need_across_edge();

      // variables required
      std::vector<Face> oneRingFaces;
      std::vector<Face> t_faces;

      // get adjacent faces
      size_t naf = this->adjacentfaces[v].size();

      if (!naf)
      {
        std::cout << "vertex " << v << " has 0 adjacent faces..." << std::endl;
      }
      else
      {
        for (int af = 0; af < naf; af++)
        {
          Face cf = this->faces[adjacentfaces[v][af]];

          t_faces.clear();
          if(IsNonObtuse(v,cf))// check angle: if non-obtuse, return existing face
          {
            t_faces.push_back(cf);
          }
          else
          {
            int nfae = this->across_edge[this->adjacentfaces[v][af]][cf.indexof(v)];
            SplitFace(t_faces,v,cf,nfae);// if obtuse, split face till we get all acute angles
          }

          for (int tf = 0; tf < t_faces.size(); tf++)
          {
            oneRingFaces.push_back(t_faces[tf]);
          }
        }
      }
      //this->colors[v] = Color::green();
      return oneRingFaces;
    }

    // FIM: initialize attributes
    //typedef std::<int> ListType;
    void InitializeAttributes(const std::vector<float>& face_speeds);

    // Debugging printout, controllable by a "verbose"ness parameter
    static int verbose;
    static void set_verbose(int);
    static int dprintf(const char *format, ...);

    // Constructor
    TriMesh() : grid_width(-1), grid_height(-1), flag_curr(0), speed_type_(ONE)
  {}
};

#endif
