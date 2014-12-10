/*
Szymon Rusinkiewicz
Princeton University

mesh_view.cc
Simple viewer
*/

#include <stdio.h>
#include <stdlib.h>
//#include "TriMesh.h"
//#include "XForm.h"
//#include "GLCamera.h"
#include <string>
#include "meshFIM.h"
#include <math.h>
#include "tetgen.h"
#include "tetmesh.h"

#include <time.h>

using std::string;

int main(int argc, char *argv[])
{
	//if (argc < 2)
	//	usage(argv[0]);

	clock_t starttime, endtime;
	
	
	  
	//const char *filename = "sphereR40_iso.ply";
	//const char *filename = "sphereR60_1k.ply";
	//const char *filename = "Atrium.ply";
	//const char *filename = "sphereR60_147237.ply";

	//const char *filename = "square.1.ply";
	//const char *filename = "dragon.ts_removedIsolated.ply";
	//const char *filename = "dragon.tx_maxSF0.5_removeIsolated.ply";
	//const char *filename = "Atrium_smoothed.ply";

    //const char *filename = "sphere_10968verts.ply";
	//const char *filename = "sphere_59021verts.ply";
	//const char *filename = "sphere_72202verts.ply";
	//const char *filename = "dragon_iso_ascii_halfsize_removeComp.ply";
	//const char *filename = "sphere_98687verts.ply";
	//const char *filename = "square.2closehole.ply";
	//const char *filename = "SquareMesh_size1024.ply";

	tetgenio in, addin, bgmin,out;
	


	//in.load_tetmesh("mouseBone.1");
	//in.load_tetmesh("long_lens260k");
	//in.load_tetmesh("lensall");
	//in.load_tetmesh("Heart_SF.5T.1_tetgen");
	//in.load_tetmesh("cube_unstruc_size256_s5");
	//in.load_tetmesh("lavalamp277k");
	in.load_tetmesh("../../levelset_data/CubeMesh_size16");
	//in.load_tetmesh("CubeMesh_size256step4");
	//int meshsize = 64;
	//in.load_tetmesh("sphereR60_1k.1");

	//const char* filename = argv[1];

	//printf("Input file name: ");
	//scanf("%s", filename);

	TetMesh themesh;
	themesh.init(in.pointlist, in.numberofpoints, in.trifacelist, in.numberoffacets,  in.tetrahedronlist, in.numberoftetrahedra, in.numberoftetrahedronattributes, in.tetrahedronattributelist );

	

	themesh.need_neighbors();
	themesh.need_adjacenttets();
	themesh.need_tet_virtual_tets();
	//themesh.need_oneringtets();
	//themesh.need_oneringstrip();


		//if (!themesh)
		//	usage(argv[0]);
		////themesh->need_normals();
		//themesh->need_tstrips();
		//themesh->need_bsphere();
		////themesh->need_across_edge();
		//meshes.push_back(themesh);

		//string xffilename = xfname(filename);
		//xffilenames.push_back(xffilename);

		//xforms.push_back(xform());
		//visible.push_back(true);
	

		meshFIM* FIMPtr = new meshFIM;

//		vector<int> seedPointList(1,0/*133152*//*20181*//*2184*/);//20181 for unstruc_s5
//		FIMPtr->SetSeedPoint(seedPointList);

		FIMPtr->SetMesh(&themesh);
		FIMPtr->FindSeedPoint();
		//FIMPtr->FindSeedPointEllipse();
		FIMPtr->InitSpeedMat();
		//FIMPtr->FindSeedPointLavalamp();

		

		int squareLength = 16;
		int squareWidth = 16;
		int squareDepth = 16;
		int squareBlockLength = 4;
		int squareBlockWidth  = 4;
		int squareBlockDepth  = 4;
		int numBlockLength = (squareLength / squareBlockLength);
		int numBlockWidth  = (squareWidth / squareBlockWidth);
		int numBlockDepth  = (squareDepth / squareBlockDepth);
		int numBlock = numBlockLength * numBlockWidth*numBlockDepth;
		int maxNumBlockVerts = squareBlockLength * squareBlockWidth * squareBlockDepth;
	
		//int maxNumBlockVerts = 64;
		//int numBlock = 64;  //10200 for longlen260k, 7500 for Heart_SF.5T.1_tetgen, 3300 for center of unstruc_s5
		FIMPtr->GraphPartition_Square(squareLength,squareWidth,squareDepth, squareBlockLength, squareBlockWidth, squareBlockDepth);
		//FIMPtr->GraphPartition_METIS2(numBlock, maxNumBlockVerts);

		FIMPtr->m_numBlock = numBlock;
		FIMPtr->PartitionTets(numBlock);
		//FIMPtr->InitializePartition(numBlock);
		FIMPtr->GenerateData();

		//endtime = clock();
		//double duration = (double)(endtime - starttime) *1000 / CLOCKS_PER_SEC;

		//printf("Computing time : %.10lf ms\n",duration);


		////FILE* vtkfile;
		////if(fopen_s(&vtkfile,"Atrium_NOISE.vtk","a+") != 0)
		////	printf_s("The vtk file was not opened\n");


		////fprintf_s(vtkfile,"POINT_DATA %d\nSCALARS traveltime float 1\nLOOKUP_TABLE default\n", FIMPtr->m_meshPtr->vertices.size());

		//float result = 0;
		//int vertIndex = 0;
		//for (int i=0;i<1/*themesh->vertices.size()*/;i++)
		//{
		//	for (int j=0;j<themesh->vertices.size();j++)
		//	{
		//		float tmp = FIMPtr->m_meshPtr->vertT[i][j];
		//		//fprintf_s(vtkfile,"%f\n",tmp);
		//		if (tmp>result)
		//		{
		//			result = tmp;
		//			vertIndex = j;
		//		}
		//		
		//	}
		//	printf("T for vert %d: %f, Vert index: %d\n", i, result, vertIndex);
		//	
		//	result = 0;
		//}

		////fclose(vtkfile);


		////float maxerror = 0.0;
		////float exactvalue = 0.0;
		////FILE* errorfile;
		////vector<float> errors;
		//////errorfile = fopen("errorfile_square16_s2.txt","w+");
		////errorfile = fopen("errorfile_SphereR60_8k_seed4511.txt","w+");

		////point seed = FIMPtr->m_meshPtr->vertices[4511];
		////point origin = point(120,120,120);
		//////point origin = point(0,0,0);
		////float Radius = 60.0;
		//////errors.resize(themesh->vertices.size());
		////float squaresum = 0;

		//////////////////coompute average radius of inscribed circle//////////

		//////FIMPtr->m_meshPtr->need_Rinscribe();
		//////
		//////vector<double> rins = FIMPtr->m_meshPtr->radiusInscribe;

		//////double rsum = 0.0;
		//////for (int i =0; i<rins.size(); i++)
		//////{
		//////	rsum+=rins[i];
		//////}

		//////double avgr = rsum / FIMPtr->m_meshPtr->faces.size();




		/////////////////////////////////////////////////////////////////////////

		////for (int i=0;i<themesh->vertices.size();i++)
		////{
		////	float computedvalue = FIMPtr->m_meshPtr->vertT[0][i];
		////	point vert = FIMPtr->m_meshPtr->vertices[i];
		////	point OS = seed - origin;
		////	point OV = vert - origin;

		////	float lenOS = sqrt( OS DOT OS);
		////	float lenOV = sqrt( OV DOT OV);

		////	float costheta = (OS DOT OV) / (lenOS*lenOV);
		////	float theta = acos(costheta);

		////	exactvalue = theta * Radius;
		////	//exactvalue = sqrt(vert[0] * vert[0] + vert[1] * vert[1]);
		////	float localerror = abs(exactvalue - computedvalue);
		////	maxerror = max(maxerror, localerror);
		////	squaresum += localerror*localerror;

		////	fprintf(errorfile,"%f\n",localerror);
		////	

		////	
		////	
		////

		////}

		////
		//////fprintf(errorfile,"avg radius of inscribed circles is: %f\n",avgr);

		////fprintf(errorfile,"max error is: %f\n",maxerror);
		////fprintf(errorfile,"rms error is: %f\n",sqrt(squaresum/(themesh->vertices.size()-1)) );
		//////fprintf(errorfile,"T for vert 0: %f, Vert index: %d\n",  result, vertIndex);
		////printf("max error is: %f\n",maxerror);
		////printf("rms error is: %f\n",sqrt(squaresum/(themesh->vertices.size()-1)));
		//////printf("avg radius of inscribed circles is: %f\n",avgr);

		////fclose(errorfile);


		////compute min max average valance 
		//int nv = FIMPtr->m_meshPtr->vertices.size();
		//int minval = 10000000, minvaltrue = 100000000;
		//int maxval = 0, maxvaltrue = 0; 
		//int sum =0, sumtrue =0;
		//for (int i =0; i<FIMPtr->m_meshPtr->vertices.size(); i++)
		//{
		//	minval = min((int)minval, (int)FIMPtr->m_meshPtr->adjacentfaces[i].size());
		//	minvaltrue = min((int)minvaltrue, (int)FIMPtr->m_meshPtr->vertOneringFaces[i].size());
		//	maxval = max((int)maxval, (int)FIMPtr->m_meshPtr->adjacentfaces[i].size());
		//	maxvaltrue = max((int)maxvaltrue, (int)FIMPtr->m_meshPtr->vertOneringFaces[i].size());
		//	sum += FIMPtr->m_meshPtr->adjacentfaces[i].size();
		//	sumtrue += FIMPtr->m_meshPtr->vertOneringFaces[i].size();


		//}

		//printf("minval is: %d\nmaxval is: %d\nminvaltrue is: %d\nmaxvaltrue is %d\naverage is %f\naveragetrue is: %f\n", minval,maxval,minvaltrue,maxvaltrue,(float)sum / (float)nv,(float)sumtrue / (float)nv);


	

		return 0;

	//glutCreateWindow("mesh_view");
	//glutDisplayFunc(redraw);
	//glutMouseFunc(mousebuttonfunc);
	//glutMotionFunc(mousemotionfunc);
	//glutKeyboardFunc(keyboardfunc);
	//glutIdleFunc(idle);

	//resetview();

	//glutMainLoop();

	


	





	
}

