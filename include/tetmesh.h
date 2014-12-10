#ifndef TETMESH_H
#define TETMESH_H
/*
Szymon Rusinkiewicz
Princeton University

TriMesh.h
Class for triangle meshes.
*/

#define  LARGENUM  10000.0
#define  SMALLNUM  0.00000001
#define  ONE       1 
#define  CURVATURE 2 
#define  NOISE     3
#define  SPEEDTYPE NOISE

#include "Vec.h"
#include "math.h"
#include <vector>
#include <list>
using std::vector;



#define MIN(a,b) ( (a)< (b) )?(a):(b)
#define MAX(a,b) ((a)>(b))?(a):(b)

class TetMesh {
//protected:
//	static bool read_helper(const char *filename, TetMesh *mesh);

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
		//float speedInv;
		float M[6];
		//float T[3];
		//Vec<3,float> edgeLens;  // edge length for 01, 12, 20

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

	//struct BBox {
	//	point min, max;
	//	point center() const { return 0.5f * (min+max); }
	//	vec size() const { return max - min; }
	//	bool valid;
	//	BBox() : valid(false)
	//		{}
	//};

	//struct BSphere {
	//	point center;
	//	float r;
	//	bool valid;
	//	BSphere() : valid(false)
	//		{}
	//};

	// Enums
	//enum tstrip_rep { TSTRIP_LENGTH, TSTRIP_TERM };
	//enum { GRID_INVALID = -1 };
	//enum speed_type { ONE = 0, CURVATURE, NOISE };

	// The basics: vertices and faces
	vector<point> vertices;
	vector<Face> faces;
	vector<Tet> tets;
	int numVert;
	int numTet;

	
	// Triangle strips
	//vector<int> tstrips;

	// Grid, if present
	//vector<int> grid;
	//int grid_width, grid_height;

	// Other per-vertex properties
	//vector<float> confidences;
	//vector<unsigned> flags;
	//unsigned flag_curr;

	// Computed per-vertex properties
	vector<vec> normals;
	vector<vec> pdir1, pdir2;
	vector<float> curv1, curv2;
	vector< Vec<4,float> > dcurv;
	vector<vec> cornerareas;
	vector<float> pointareas;
	vector< vector<float> > vertT;

	// Bounding structures
	//BBox bbox;
	//BSphere bsphere;

	// Connectivity structures:
	//  For each vertex, all neighboring vertices
	vector< vector<int> > neighbors;
	//  For each vertex, all neighboring faces
	vector< vector<int> > adjacenttets;
	vector< vector<int> > oneringstrips;
	vector< vector<float> > oneringspeedI;

	vector<double> radiusInscribe;
	vector<Tet> across_face;
	


	vector< vector<Tet> > vertOneringTets;
	vector< vector<Tet> > tetVirtualTets;
	//  For each face, the three faces attached to its edges
	//  (for example, across_edge[3][2] is the number of the face
	//   that's touching the edge opposite vertex 2 of face 3)
	//vector<Face> across_edge;

	vector<float> noiseOnVert;
	//vector<float> noiseOnFace;
	

	//int SPEEDTYPE;
	// Compute all this stuff...
	//void setSpeedType(int st)
	//{
		//ST = st;
	//}
	//void need_tstrips();
	//void convert_strips(tstrip_rep rep);
	//void unpack_tstrips();
	//void triangulate_grid();
	//void need_faces()
	//{
	//	if (!faces.empty())
	//		return;
	//	if (!tstrips.empty())
	//		unpack_tstrips();
	//	else if (!grid.empty())
	//		triangulate_grid();
	//}

	void need_faceedges();
	void need_speed();
	void need_noise();
	void need_oneringtets();
	void need_normals();
	void need_pointareas();
	//void need_curvatures();
	//void need_dcurv();
	//void need_bbox();
	//void need_bsphere();
	void need_neighbors();
	void need_adjacenttets();
	void need_oneringstrip();
	//void need_across_edge();
	void need_meshinfo();
	void need_Rinscribe();
	void need_across_face();

	bool IsNonObtuse(int v, Tet t);
	void SplitFace(vector<Tet> &acTets, int v, Tet ct, int nfAdj);
	vector<Tet> GetOneRing(int v);
	void need_tet_virtual_tets();

	void InitializeAttributes(int currentVert , vector<int> seeds )
	{
		// initialize the travel times over all vertices...
		int nv = this->vertices.size();

		for (int v = 0; v < nv; v++)
		{			
			this->vertT[currentVert].push_back(LARGENUM);  //modified from this->vertT[v] = 1000000.0)
		}

		//vector<int> nb;

		// initialize seed points if present...
		if (!seeds.empty())
		{
			int ns = seeds.size();
			for (int s = 0; s < ns; s++)
			{
				this->vertT[currentVert][seeds[s]] = 0.0;  //modified from this->vertT[s] = 0.0;
				//nb = this->neighbors[seeds[s]];
				
			}


		}
	}

	// Debugging printout, controllable by a "verbose"ness parameter
	static int verbose;
	static void set_verbose(int);
	static int dprintf(const char *format, ...);

	
	void init(double* pointlist, int numpoint, int*trilist, int numtri, int* tetlist, int numtet, int numattr, double* attrlist);

	//Constructor
	TetMesh(){}
};

#endif
