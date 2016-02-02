#ifndef TETMESH_H
#define TETMESH_H
/*
   Szymon Rusinkiewicz
   Princeton University

   TriMesh.h
   Class for triangle meshes.
 */

#ifndef  LARGENUM
#define  LARGENUM  1000000000
#endif
#ifndef  SMALLNUM
#define  SMALLNUM  0.00000001
#endif

#ifndef ONE
#define ONE 1
#endif
#ifndef NOISE
#define NOISE 2
#endif

#include "Vec.h"
#include "math.h"
#include <vector>
#include <list>

#ifndef M_PI
#define M_PI 3.14159265359
#endif

class TetMesh {
  //protected:
  //  static bool read_helper(const char *filename, TetMesh *mesh);

  public:
    // Types
    struct Face{
      int v[3];

      float speedInv;
      float T[3];
      Vec<3,float> edgeLens;  // edge length for 01, 12, 20

      Face() {}
      Face(const int &v0, const int &v1, const int &v2)
      {
        v[0] = v0; v[1] = v1; v[2] = v2;
      }
      Face(const int *v_)
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



    struct Tet {
      int v[4];
      int obtuseV;
      float M[6];

      Tet() {}
      Tet(const int &v0, const int &v1, const int &v2, const int &v3)
      {

        v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;

      }
      Tet(const int *v_)
      {
        v[0] = v_[0]; v[1] = v_[1]; v[2] = v_[2]; v[3] = v_[3];

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
          (v[2] == v_) ? 2 :
          (v[3] == v_) ? 3 : -1;
      }
    };

    // The basics: vertices and faces
    std::vector<point> vertices;
    std::vector<Face> faces;
    std::vector<Tet> tets;
    int numVert;
    int numTet;


    // Computed per-vertex properties
    std::vector<vec> normals;
    std::vector<vec> pdir1, pdir2;
    std::vector<float> curv1, curv2;
    std::vector< Vec<4,float> > dcurv;
    std::vector<vec> cornerareas;
    std::vector<float> pointareas;
    std::vector<float> vertT;

    // Connectivity structures:
    //  For each vertex, all neighboring vertices
    std::vector< std::vector<int> > neighbors;
    //  For each vertex, all neighboring faces
    std::vector< std::vector<int> > adjacenttets;
    std::vector< std::vector<int> > oneringstrips;
    std::vector< std::vector<float> > oneringspeedI;

    std::vector<float> radiusInscribe;
    std::vector<Tet> across_face;



    std::vector< std::vector<Tet> > vertOneringTets;
    std::vector< std::vector<Tet> > tetVirtualTets;

    std::vector<float> noiseOnVert;

    void reorient();
    void need_faceedges();
    void need_speed(int speedtype);
    void need_noise(float low, float high);
    void need_oneringtets();
    void need_normals();
    void need_pointareas();
    void need_neighbors(bool verbose = false);
    void need_adjacenttets(bool verbose = false);
    void need_oneringstrip();
    void need_meshinfo();
    void need_Rinscribe();
    void need_across_face(bool verbose = false);

    bool IsNonObtuse(int v, Tet t);
    void SplitFace(std::vector<Tet> &acTets, int v, Tet ct, int nfAdj);
    std::vector<Tet> GetOneRing(int v);
    void need_tet_virtual_tets(bool verbose = false);

    // Debugging printout, controllable by a "verbose"ness parameter
    static int verbose;
    static void set_verbose(int);
    static int dprintf(const char *format, ...);


    void init(float* pointlist, int numpoint, int*trilist, int numtri,
      int* tetlist, int numtet, float* attrlist, 
      std::vector<float> speedMtx, 
      bool verbose = false);

    //Constructor
    TetMesh(){}
};

#endif
