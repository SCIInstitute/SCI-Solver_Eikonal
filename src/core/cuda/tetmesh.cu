#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "tetmesh.h"
//#include "Stdafx.h"

using namespace std;


void TetMesh::init(double* pointlist, int numpoint, int*trilist, int numtri, int* tetlist, int numtet, int numattr, double* attrlist)
{
	numVert = numpoint;
	numTet  = numtet;

	printf("Total num of vertices is %d\n", numVert);
	printf("Total num of tets is     %d\n", numTet);
	vertices.resize(numpoint);
	tets.resize(numtet);
	for(int i =0; i< numpoint; i++)
	{
		vertices[i][0] = pointlist[3*i+0];
		vertices[i][1] = pointlist[3*i+1];
		vertices[i][2] = pointlist[3*i+2];




	}

	//check the index start from 0 or 1
	int minidx = 1000000000;
	for(int i=0; i<numtet*4; i++)
	{
		minidx = MIN(minidx, tetlist[i]);

	}

	if(minidx == 0)
	{
		for(int i =0; i< numtet; i++)
		{

			tets[i][0] = tetlist[4*i+0];  
			tets[i][1] = tetlist[4*i+1];
			tets[i][2] = tetlist[4*i+2];
			tets[i][3] = tetlist[4*i+3];
			tets[i].obtuseV = -1;
			tets[i].M[0] = 1.0;
			tets[i].M[1] = 0.0;
			tets[i].M[2] = 0.0;
			tets[i].M[3] = 1.0;
			tets[i].M[4] = 0.0;
			tets[i].M[5] = 1.0;

		}
	}
	else if(minidx == 1)
	{
		for(int i =0; i< numtet; i++)
		{
			tets[i][0] = tetlist[4*i+0]-1;  // -1 because the oringal index is from 1 and change it to 0
			tets[i][1] = tetlist[4*i+1]-1;
			tets[i][2] = tetlist[4*i+2]-1;
			tets[i][3] = tetlist[4*i+3]-1;
			tets[i].obtuseV = -1;
			tets[i].M[0] = 1.0;
			tets[i].M[1] = 0.0;
			tets[i].M[2] = 0.0;
			tets[i].M[3] = 1.0;
			tets[i].M[4] = 0.0;
			tets[i].M[5] = 1.0;
		}
	}
	else
	{
		printf("error!!! index not start from 0 or 1!!\n");
	}


	if(numattr > 0)
	{
		for(int i =0; i< numtet; i++)
		{
			int mat = (int)attrlist[i];
			switch(mat)
			{
			case 0:
				tets[i].M[0] = 1.0;
				tets[i].M[1] = 0.0;
				tets[i].M[2] = 0.0;
				tets[i].M[3] = 1.0;
				tets[i].M[4] = 0.0;
				tets[i].M[5] = 1.0;
				break;
			case 1:
				tets[i].M[0] = 10.0;
				tets[i].M[1] = 0.0;
				tets[i].M[2] = 0.0;
				tets[i].M[3] = 10.0;
				tets[i].M[4] = 0.0;
				tets[i].M[5] = 10.0;
				break;


			}
		}


	}
	else
	{
		for(int i =0; i< numtet; i++)
		{
			tets[i].M[0] = 1.0;
			tets[i].M[1] = 0.0;
			tets[i].M[2] = 0.0;
			tets[i].M[3] = 1.0;
			tets[i].M[4] = 0.0;
			tets[i].M[5] = 1.0;
		}

	}


	//if(numattr > 0)
	//{
	//	for(int i =0; i< numtet; i++)
	//	{
	//		int mat = (int)attrlist[i];
	//		switch(mat)
	//		{
	//		case 1:
	//			tets[i].speedInv = 1.0;
	//			break;
	//		case 2:
	//			tets[i].speedInv = 1.52; // refractive index of glass 1.52
	//			break;
	//		case 3:
	//			tets[i].speedInv = 1.0;
	//			break;


	//		}
	//	}


	//}
	//else
	//{
	//	for(int i =0; i< numtet; i++)
	//		tets[i].speedInv = 1.0;

	//}


	//for(int i =0; i< numtet; i++)
	//{
	//	tets[i].f[0] = Face(tets[i][0], tets[i][1], tets[i][2]);
	//	tets[i].f[1] = Face(tets[i][1], tets[i][2], tets[i][3]);
	//	tets[i].f[2] = Face(tets[i][2], tets[i][3], tets[i][0]);
	//	tets[i].f[3] = Face(tets[i][3], tets[i][0], tets[i][1]);





	//}


}
// Find the direct neighbors of each vertex
void TetMesh::need_neighbors()
{
	if (!neighbors.empty())
		return;


	cout << "Finding vertex neighbors... " << endl;
	int nv = vertices.size(), nt = tets.size();

	//vector<int> numneighbors(nv);
	//for (int i = 0; i < nt; i++) {
	//	numneighbors[tets[i][0]]++;
	//	numneighbors[tets[i][1]]++;
	//	numneighbors[tets[i][2]]++;
	//	numneighbors[tets[i][3]]++;
	//}

	neighbors.resize(nv);
	//for (int i = 0; i < nv; i++)
	//	neighbors[i].reserve(numneighbors[i]+2); // Slop for boundaries

	for (int i = 0; i < nt; i++) {
		for (int j = 0; j < 4; j++) {
			vector<int> &me = neighbors[tets[i][j]];
			int n1 = tets[i][(j+1)%4];
			int n2 = tets[i][(j+2)%4];
			int n3 = tets[i][(j+3)%4];
			if (find(me.begin(), me.end(), n1) == me.end())
				me.push_back(n1);
			if (find(me.begin(), me.end(), n2) == me.end())
				me.push_back(n2);
			if (find(me.begin(), me.end(), n3) == me.end())
				me.push_back(n3);
		}
	}

	cout << "Done.\n" << endl;
}


// Find the tets touching each vertex
void TetMesh::need_adjacenttets()
{
	if (!adjacenttets.empty())
		return;

	cout << "Finding vertex to triangle maps... " << endl;
	int nv = vertices.size(), nt = tets.size();

	//vector<int> numadjacentfaces(nv);
	//for (int i = 0; i < nt; i++) {
	//	numadjacentfaces[tets[i][0]]++;
	//	numadjacentfaces[tets[i][1]]++;
	//	numadjacentfaces[tets[i][2]]++;
	//}

	adjacenttets.resize(vertices.size());
	//for (int i = 0; i < nv; i++)
	//	adjacentfaces[i].reserve(numadjacentfaces[i]);

	for (int i = 0; i < nt; i++) {
		for (int j = 0; j < 4; j++)
			adjacenttets[tets[i][j]].push_back(i);
	}

	cout << "Done.\n" << endl;
}


void TetMesh::need_speed()
{
}


//void TetMesh::need_oneringstrip()
//{
//	
//
//	int numVert = vertices.size();
//	//oneringstrip.clear();
//	oneringstrips.resize(numVert);
//	oneringspeedI.resize(numVert);
//	for(int i =0; i<numVert; i++)
//	{
//		udword* topology;
//		vector<int> adjtets = adjacenttets[i];
//		//vector<Tet> adjtets = vertOneringTets[i];
//		topology = (udword*)malloc(3*sizeof(udword)*adjtets.size());
//		for(int j=0; j<adjtets.size(); j++)
//		{
//			int tmp=0;
//			for(int k=0; k<4; k++)
//			{
//				
//				if(tets[adjtets[j]][k] != i)
//				//if(adjtets[j][k] != i)
//				{
//					topology[j*3+tmp] = tets[adjtets[j]][k];
//					//topology[j*3+tmp] = adjtets[j][k];
//					tmp++;
//				}
//
//
//			}
//
//		}
//
//		STRIPERCREATE sc;
//		sc.DFaces			= topology;
//		sc.NbFaces			= adjtets.size();
//		sc.AskForWords		= false;
//		sc.ConnectAllStrips	= true;
//		sc.OneSided			= false;
//		sc.SGIAlgorithm		= true;
//
//		Striper Strip;
//		Strip.Init(sc);
//
//		STRIPERRESULT sr;
//		Strip.Compute(sr);
//
//		vector< vector<int> > strips;
//		strips.resize(sr.NbStrips);
//		udword* Refs = (udword*)sr.StripRuns;
//
//		oneringstrips[i].resize( sr.NbStrips);
//		
//		for(udword k=0;k<sr.NbStrips;k++)
//		{
//			
//			udword NbRefs = sr.StripLengths[k];
//			for(udword j=0;j<NbRefs;j++)
//			{
//				strips[k].push_back(*Refs++);
//				//oneringstrip[i][k].push_back(*Refs++);
//			}
//			
//		}
//
//		
//		//printf("vert %d num strips: %d\n",i, sr.NbStrips);
//
//		
//		oneringstrips[i] = strips[0];
//		
//
//
//
//
//
//		free(topology);
//
//
//
//
//
//	}
//
//	
//
//}

bool TetMesh::IsNonObtuse(int v, Tet t)
{
	int D = t.indexof(v);
	int A = (D+1)%4; 
	int B = (D+2)%4; 
	int C = (D+3)%4; 

	point P1 = vertices[t[A]];
	point P2 = vertices[t[B]];
	point P3 = vertices[t[C]];
	point P4 = vertices[t[D]];

	point a = P1 - P4;
	point b = P2 - P4;
	point c = P3 - P4;

	float det = fabs((a ^ (b % c)));

	float al = a.norm();
	float bl = b.norm();
	float cl = c.norm();

	float div = al*bl*cl + (a ^ b)*cl + (a ^ c)*bl + (b ^ c)*al;
	float at = atan2(det, div);
	if(at < 0) at += M_PI; // If det>0 && div<0 atan2 returns < 0, so add pi.

	//return omega < M_PI / 2.0;
	return 1;

}


void TetMesh::SplitFace(vector<Tet> &acTets, int v, Tet ct, int nfAdj)
{		
	// get all the four vertices in order
	/* v1         v4
	+-------+
	\     . \
	\   .   \
	\ .     \
	+-------+
	v2         v3 */

	need_neighbors();
	int iV = ct.indexof(v);											// get index of v in terms of cf
	int v1 = v;
	int v2 = ct[(iV+1)%4];
	int v3 = ct[(iV+2)%4];
	int v4 = ct[(iV+3)%4];
	iV = tets[nfAdj].indexof(v2);				// get index of v in terms of adjacent face

	int v5;
	for(int i=0; i<4; i++)
	{
		if(tets[nfAdj][i] != v2 && tets[nfAdj][i] != v3 && tets[nfAdj][i] != v4)
			v5 = tets[nfAdj][i];

	}
	neighbors[v5].push_back(v1);
	//Tet af = tets[nfAdj];

	// create faces (v1,v3,v4) and (v1,v2,v3), check angle at v1
	Tet t1(v1, v2, v3, v5);	
	Tet t2(v1, v3, v4, v5);	
	Tet t3(v1, v2, v4, v5);

	

	if (IsNonObtuse(v,t1))
	{
		acTets.push_back(t1);
	}
	else
	{
		int nfAdj_new = across_face[nfAdj][tets[nfAdj].indexof(v4)];
		if (nfAdj_new > -1)
		{
			SplitFace(acTets,v,t1,nfAdj_new);
			
		}
		else
			printf("NO cross edge!!! Maybe a hole!!\n");
		//SplitFace(acFaces,v,f1,nfAdj_new, currentVert);
	}

	if (IsNonObtuse(v,t2))
	{
		acTets.push_back(t2);
	}
	else
	{
		int nfAdj_new = across_face[nfAdj][tets[nfAdj].indexof(v2)];
		if (nfAdj_new > -1)
		{
			SplitFace(acTets,v,t2,nfAdj_new/*,currentVert*/);
		}
		else
			printf("NO cross edge!!! Maybe a hole!!\n");
		//SplitFace(acFaces,v,f2,nfAdj_new,currentVert);
	}

	if (IsNonObtuse(v,t3))
	{
		acTets.push_back(t3);
	}
	else
	{
		int nfAdj_new = across_face[nfAdj][tets[nfAdj].indexof(v3)];
		if (nfAdj_new > -1)
		{
			SplitFace(acTets,v,t3,nfAdj_new/*,currentVert*/);
		}
		else
			printf("NO cross edge!!! Maybe a hole!!\n");
		//SplitFace(acFaces,v,f2,nfAdj_new,currentVert);
	}
}


void TetMesh::need_across_face()
{
	if (!across_face.empty())
		return;
	need_adjacenttets();

	printf("Finding across-face maps... ");

	int nt = tets.size();
	across_face.resize(nt, Tet(-1,-1,-1, -1));

	for (int i = 0; i < nt; i++) {
		for (int j = 0; j < 4; j++) {
			if (across_face[i][j] != -1)
				continue;
			int v1 = tets[i][(j+1)%4];
			int v2 = tets[i][(j+2)%4];
			int v3 = tets[i][(j+3)%4];
			const vector<int> &a1 = adjacenttets[v1];
			const vector<int> &a2 = adjacenttets[v2];
			const vector<int> &a3 = adjacenttets[v3];
			for (int k1 = 0; k1 < a1.size(); k1++) 
			{
				int other = a1[k1];
				if (other == i)
					continue;
				vector<int>::const_iterator it =
					find(a2.begin(), a2.end(), other);

				vector<int>::const_iterator it2 =
					find(a3.begin(), a3.end(), other);

				if (it == a2.end() || it2 == a3.end())
					continue;

				across_face[i][j] = other;
				break;

				
				
			}
		}
	}

	printf("Done.\n");
}

vector<TetMesh::Tet> TetMesh::GetOneRing(int v)
{
	// make sure we have the across-edge map
	if (across_face.empty())
		need_across_face();

	// variables required
	vector<Tet> oneRingTets;
	vector<Tet> t_tets;

	// get adjacent faces
	int naf = adjacenttets[v].size();

	if (!naf)
	{
		std::cout << "vertex " << v << " has 0 adjacent faces..." << std::endl;
	}
	else
	{
		for (int af = 0; af < naf; af++)
		{
			Tet ct = this->tets[adjacenttets[v][af]];

			t_tets.clear();
			if(IsNonObtuse(v,ct))// check angle: if non-obtuse, return existing face
			{
				//this->colors[cf[0]] = Color::red();
				//this->colors[cf[1]] = Color::red();
				//this->colors[cf[2]] = Color::red();
				//t_tets.push_back(ct);					
			}
			else
			{
				int nfae = this->across_face[adjacenttets[v][af]][ct.indexof(v)];
				if (nfae > -1)
				{
					SplitFace(t_tets,v,ct,nfae/*,currentVert*/);// if obtuse, split face till we get all acute angles
				}
				else
					printf("NO cross edge!!! Maybe a hole!!\n");
				//SplitFace(t_faces,v,cf,nfae,currentVert);// if obtuse, split face till we get all acute angles
			}

			for (int tf = 0; tf < t_tets.size(); tf++)
			{					
				//this->colors[t_faces[tf][0]] = Color::red();
				//this->colors[t_faces[tf][1]] = Color::red();
				//this->colors[t_faces[tf][2]] = Color::red();
				oneRingTets.push_back(t_tets[tf]);
			}
		}				
	}
	//this->colors[v] = Color::green();
	return oneRingTets;
}

void TetMesh::need_oneringtets()
{

	if (vertOneringTets.empty())
	{
		vertOneringTets.resize(vertices.size());
		for (int i=0; i< vertices.size();i++)
		{
			vertOneringTets[i] = GetOneRing(i);

		}
		
	}
	
}

void TetMesh::need_tet_virtual_tets()
{

	
	need_across_face();
	vector<Tet> t_tets;
	Tet t;
	int numTets = tets.size();
	tetVirtualTets.resize(numTets);
	for (int i = 0; i < numTets; i++)
	{
		t_tets.clear();
		t = tets[i];

		for (int j = 0; j< 4 ; j++)
		{
			if(!IsNonObtuse(t[j],t))// check angle: if non-obtuse, return existing face
			{
				t.obtuseV = j;

				int nfae = across_face[i][j];
				if (nfae > -1)
				{
					SplitFace(t_tets,t[j],t,nfae);// if obtuse, split face till we get all acute angles
				}
				else
					printf("NO cross edge!!! Maybe a hole!!\n");
				
			}
		}

		tetVirtualTets[i] = t_tets;
	}
	

	
	
}
