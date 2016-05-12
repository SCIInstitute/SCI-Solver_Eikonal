GPUTUM: Unstructured Eikonal
=====================

GPUTUM: Unstructured Eikonal is a C++/CUDA library written to solve the Eikonal equation on a 3D tet meshes and 2D triangular meshes. It uses the fast iterative method (FIM) to solve efficiently, and uses GPU hardware.

The code was written by Zhisong Fu at the Scientific Computing and Imaging Institute, 
University of Utah, Salt Lake City, USA. The theory behind this code is published in the papers linked below. 
Table of Contents
========
- [Eikonal 2D Aknowledgements](#eikonal-2d-aknowledgements)
- [Eikonal 3D Aknowledgements](#eikonal-3d-aknowledgements)
- [Requirements](#requirements)
- [Building](#building)<br/>
		- [Linux / OSX](#linux-and-osx)<br/>
		- [Windows](#windows)<br/>
- [Running Examples](#running-examples)
- [Using the Library](#using-the-library)
- [Testing](#testing)<br/>

<h4>Eikonal 2D Aknowledgements</h4>
**<a href ="http://epubs.siam.org/doi/abs/10.1137/100788951">A Fast Iterative Method for Solving the 
Eikonal Equation on Triangulated Surfaces</a>**<br/>
<img src="https://raw.githubusercontent.com/SCIInstitute/SCI-Solver_Eikonal/master/src/Resources/eikonal2d.png"  align="right" hspace="20" width=450>
*NOTE*<br/>
This code is a merge of Zhisong's 2D and 3D Unstructured Eikonal code. The oldest commit history for 2D can be found <a href="https://github.com/SCIInstitute/SCI-Solver_Eikonal-2D">here</a>.

**AUTHORS:**
<br/>Zhisong Fu(*a*) <br/>
Won-Ki Jeong(*b*) <br/>
Yongsheng Pan(*a*) <br/>
Robert M. Kirby(*a*) <br/>
Ross T. Whitaker(*a*) <br/>

This library solves for the Eikional values on vertices located on a triangular mesh. Several mesh formats
are supported, and are read by the <a href="http://graphics.stanford.edu/software/trimesh/">TriMesh Library</a>. 
The <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/download">METIS library</a> is used to partition unstructured 
meshes. <a href="https://code.google.com/p/googletest/">
Google Test</a> is used for testing.
<br/>

<h4>Eikonal 3D Aknowledgements</h4>
**<a href="http://epubs.siam.org/doi/abs/10.1137/120881956"> A Fast Iterative Method for 
Solving the Eikonal Equation on Tetrahedral Domains</a>**<br/>
<img src="https://raw.githubusercontent.com/SCIInstitute/SCI-Solver_Eikonal/master/src/Resources/eikonal3d.png"  align="right" hspace="20" width=450>

**AUTHORS:**
Zhisong Fu(a,b) <br/>
Robert M. Kirby(a,b) <br/>
Ross T. Whitaker(a,b) <br/>

This library solves for the Eikional values on vertices located on a tetrahedral mesh. Several mesh formats
are supported, and are read by the <a href="http://wias-berlin.de/software/tetgen/">TetGen Library</a>. 
The <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/download">METIS library</a> is used to partition unstructured 
meshes. <a href="https://code.google.com/p/googletest/">
Google Test</a> is used for testing.
<br/><br/><br/><br/>
Requirements
==============

 * Git, CMake (2.8+ recommended), and the standard system build environment tools.
 * You will need a CUDA Compatible Graphics card. See <a href="https://developer.nvidia.com/cuda-gpus">here</a> You will also need to be sure your card has CUDA compute capability of at least 2.0.
 * SCI-Solver_Eikonal is compatible with the latest CUDA toolkit (7.5). Download <a href="https://developer.nvidia.com/cuda-downloads">here</a>.
 * This project has been tested on OpenSuse 13.1 (Bottle) on NVidia GeForce GTX 680 HD, Windows 7 on NVidia GeForce GTX 775M, and OSX 10.10 on NVidia GeForce GTX 775M. 
 * If you have a CUDA graphics card equal to or greater than our test machines and are experiencing issues, please contact the repository owners.
 * Windows: You will need Microsoft Visual Studio 2010+ build tools. This document describes the "NMake" process.
 * OSX: Please be sure to follow setup for CUDA <a href="http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-mac-os-x/#axzz3W4nXNNin">here</a>. There are several compatability requirements for different MAC machines, including using a different version of CUDA (ie. 5.5).

Building
==============

<h3>Linux and OSX</h3>
In a terminal:
```c++
mkdir SCI-SOLVER_Eikonal/build
cd SCI-SOLVER_Eikonal/build
cmake ../src
make
```

<h3>Windows</h3>
Open a Visual Studio (32 or 64 bit) Native Tools Command Prompt. 
Follow these commands:
```c++
mkdir C:\Path\To\SCI-Solver_Eikonal\build
cd C:\Path\To\SCI-Solver_Eikonal\build
cmake -G "NMake Makefiles" ..\src
nmake
```

**Note:** For all platforms, you may need to specify your CUDA toolkit location (especially if you have multiple CUDA versions installed):
```c++
cmake -DCUDA_TOOLKIT_ROOT_DIR="~/NVIDIA/CUDA-7.5" ../src
```
(Assuming this is the location).

**Note:** If you have compile errors such as <code>undefined reference: atomicAdd</code>, it is likely you need to set your compute capability manually. CMake outputs whether compute capability was determined automatically, or if you need to set it manually. The default minimum compute capability is 2.0.

```c++
cmake -DCUDA_COMPUTE_CAPABILITY=20 ../src
make
```

Running Examples
==============

You will need to enable examples in your build to compile and run them.

```c++
cmake -DBUILD_EXAMPLES=ON ../src
make
```

You will find the example binaries built in the <code>build/examples</code> directory.

Run the examples in the build directory:

```c++
examples/Example1 
examples/Example2  
...
```
Each example has a <code>-h</code> flag that prints options for that example. <br/>

Follow the example source code in <code>src/examples</code> to learn how to use the library.
<br/>
To run examples similar to the paper, the following example calls would do so:<br/>
<b>2D Eikonal Homogeneous, Sphere </b>(add <code>-b 10</code> or another #verts/block if 2D fails)<br/>
<code>examples/Example2 -v -i ../src/test/test_data/sphere_1154verts.ply</code><br/>
<b>2D Eikonal Scalar Heterogeneous, Square </b> <br/>
<code>examples/Example2 -v -i ../src/test/test_data/SquareMesh_size32.ply -x ../src/test/test_data/square_scalars.txt</code><br/>
<b>3D Eikonal Homogeneous, Sphere </b><br/>
<code>examples/Example1 -v -i ../src/test/test_data/sphere8092</code><br/>
<b>3D Eikonal Tensor Heterogeneous, Torus </b> <br/>
<code>examples/Example1 -v -i ../src/test/test_data/torus -x ../src/test/test_data/torus_matrices.txt -n 300</code><br/>
<b>3D Eikonal Scalar Heterogeneous, Cube </b> <br/>
<code>examples/Example1 -v -i ../src/test/test_data/CubeMesh_size256step8_correct -x ../src/test/test_data/cube_scalars.txt -n 300</code><br/>

**NOTE** All examples output a <code>result0.vtk</code> VTK file in the current directory.
These files are easily viewed via VTK readers like Paraview. You can clip and add iso-values to 
more distinctly visualize the result.

Using the Library
==============

A basic usage of the library links to the <code>EIKONAL_CORE</code> 
library during build and includes the headers needed, which are usually no more than:

```c++
#include <Eikonal.h>
```

Then a program would setup the Eikonal parameters using the 
<code>Eikonal object</code> object and call 
<code>object.solveEikonal()</code> to generate
the array of vertex values per iteration.

Here is a minimal usage example (in 3D).<br/>
```c++
#include <Eikonal.h>
#include <iostream>
int main(int argc, char *argv[])
{
  Eikonal data(false);
  //the below means ~/my_tet_mesh.node & ~/my_tet_mesh.ele
  data.filename_ = "~/my_tet_mesh"; 
  //Run the solver
  data.solveEikonal();
  //now use the result
  for(size_t i = 0; i < data.getFinalResult().size(); i++)
    std::cout << "Vertex " << i << " value: " << 
      data.getFinalResult()[i] << std::endl;
  return 0;
}
```

The following helper functions are available before running the solver:
```c++
void Eikonal::initializeMesh();
void Eikonal::initSpeedMtxMultipliers(std::vector<float> speeds); //set the speeds per face/tet
void Eikonal::initializeVertices(std::vector<float> values); //set the seed points to 0, others to LARGENUM
```
The following helper functions are available after running the solver:
```c++
std::vector < float > Eikonal::getFinalResult();
std::vector < float > Eikonal::getResultAtIteration(size_t i);
size_t Eikonal::numIterations(); 
void Eikonal::writeVTK(bool writeEachIteration); //write the vtk files that show progression of the solver
```
You can also access the results and the mesh directly after running the solver:
```c++
TetMesh * Eikonal::tetMesh_;
// OR
TriMesh * Eikonal::triMesh_;
// AND
std::vector < std::vector < float > > Eikonal::iteration_values_;
```

<h3>Eikonal Options</h3>

```C++
  class Eikonal {
      bool verbose_;                    //option to set for runtime verbosity [Default false]
      std::string filename_;            //the input tet mesh filename         [Default ../src/test/test_data/sphere339
      int maxBlocks_;                   //the max # of blocks (patches)
                                        //   on the convergence queue         [Default 100]
      int maxVertsPerBlock_;            //Max # of vertices per block         [Default 64]
      float stopDistance_;            	//Stop distance for Eikonal	      [Default 50000.]
      bool isStructured_;               //Whether the mesh is structured      [Default false]
      int squareLength_;                //if structured, the square size      [Default 16, 16, 16]
      int squareWidth_;                 
      int squareDepth_;                 
      int squareBlockLength_;           //if structured, CUDA block length    [Default 1]
      int squareBlockWidth_;            //if structured, CUDA block width     [Default 1]
      int squareBlockDepth;             //if structured, CUDA block width     [Default 1]
      int maxIterations_;               //when to stop iterating if fail      [Default 100]
      bool isTriMesh_;                  //This is a triangle mesh             [Default true]
  };
```
<br/>
You will need to make sure your CMake/Makfile/Build setup knows where 
to point for the library and header files. See the examples and their CMakeLists.txt.<br/><br/>

**NOTE** The maxVertsPerBlock_ parameter for the 2D library often needs to be tuned per mesh. If the library fails
to converge all vertices (reaches the max iterations), try a different number for this parameter. Various 
possiblities per mesh ranging from 2-2000 may improve error once a good range is determined.<br/>

Testing
==============
The repo comes with a set of regression tests to see if recent changes break 
expected results. To build the tests, you will need to set 
<code>BUILD_TESTING</code> to "ON" in either <code>ccmake</code> or when calling CMake:

```c++
cmake -DBUILD_TESTING=ON ../src
```
After building, run <code>make test</code> or <code>ctest</code> in the build directory to run tests.<br/>
<h4>Windows</h4>
The gtest library included in the repo needs to be built with 
forced shared libraries on Windows, so use the following:

```c++
cmake -DBUILD_TESTING=ON -Dgtest_forced_shared_crt=ON ../src
```
Be sure to include all other necessary CMake definitions as annotated above.
