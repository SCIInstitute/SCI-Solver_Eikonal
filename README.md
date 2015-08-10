SCI-Solver_Eikonal
=====================

SCI-Solver_Eikonal is a C++/CUDA library written to solve the Eikonal equation on a 3D tet meshes and 2D triangular meshes. It uses the fast iterative method (FIM) to solve efficiently, and uses GPU hardware.

The code was written by Zhisong Fu at the Scientific Computing and Imaging Institute, 
University of Utah, Salt Lake City, USA. The theory behind this code is published in the papers: 

**<a href ="http://epubs.siam.org/doi/abs/10.1137/100788951">A Fast Iterative Method for Solving the 
Eikonal Equation on Triangulated Surfaces</a>**
<img src="https://raw.githubusercontent.com/SCIInstitute/SCI-Solver_Eikonal/master/src/Resources/eikonal2d.png"  align="right" hspace="20">

**AUTHORS:** <br/>
<br/>Zhisong Fu(*a*) <br/>
Won-Ki Jeong(*b*) <br/>
Yongsheng Pan(*a*) <br/>
Robert M. Kirby(*a*) <br/>
Ross T. Whitaker(*a*) <br/>

This library solves for the Eikional values on vertices located on a triangular mesh. Several mesh formats
are supported, and are read by the <a href="http://graphics.stanford.edu/software/trimesh/">TriMesh Library</a>. 
The <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/download">METIS library</a> is used to partition unstructured 
meshes. Otherwise, limited structured meshes are supported. <a href="https://code.google.com/p/googletest/">
Google Test</a> is used for testing.

**<a href="http://epubs.siam.org/doi/abs/10.1137/120881956"> A Fast Iterative Method for 
Solving the Eikonal Equation on Tetrahedral Domains</a>**

**AUTHORS:** <br/>
Zhisong Fu(a,b) <br/>
Robert M. Kirby(a,b) <br/>
Ross T. Whitaker(a,b) <br/>

Requirements
==============

 * Git, CMake (3.0+ recommended), and the standard system build environment tools.
 * You will need a CUDA Compatible Graphics card. See <a href="https://developer.nvidia.com/cuda-gpus">here</a> You will also need to be sure your card has CUDA compute capability of at least 2.0.
 * SCI-Solver_Eikonal is compatible with the latest CUDA toolkit (7.0). Download <a href="https://developer.nvidia.com/cuda-downloads">here</a>.
 * This project has been tested on OpenSuse 13.1 (Bottle) on NVidia GeForce GTX 680 HD, Windows 7 on NVidia GeForce GTX 775M, and OSX 10.10 on NVidia GeForce GTX 775M. 
 * If you have a CUDA graphics card equal to or greater than our test machines and are experiencing issues, please contact the repository owners.
 * Windows: You will need Microsoft Visual Studio 2010+ build tools. This document describes the "NMake" process.
 * OSX: Please be sure to follow setup for CUDA <a href="http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-mac-os-x/#axzz3W4nXNNin">here</a>. There are several compatability requirements for different MAC machines, including using a different version of CUDA (ie. 5.5).

Building
==============

<h3>Unix / OSX</h3>
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
cmake -DCUDA_TOOLKIT_ROOT_DIR="~/NVIDIA/CUDA-7.0" ../src
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

Using the Library
==============

A basic usage of the library links to the <code>Eikonal-2D_CORE</code> or  <code>Eikonal-3D_CORE</code> 
library during build and includes the headers needed, which are usually no more than:

```c++
#include <Eikonal2D.h>
// OR
#include <Eikonal3D.h>
```

Then a program would setup the Eikonal parameters using the 
<code>Eikonal::Eikonal2D -OR- Eikonal::Eikonal3D</code> object and call 
<code>Eikonal::solveEikonal2D() -OR- Eikonal::solveEikonal3D()</code> to generate
the array of vertex values per iteration.

You will need to make sure your CMake/Makfile/Build setup knows where 
to point for the library and header files. See the examples and their CMakeLists.txt.

Testing
==============
The repo comes with a set of regression tests to see if recent changes break 
expected results. To build the tests, you will need to set 
<code>BUILD_TESTING</code> to "ON" in either <code>ccmake</code> or when calling CMake:

```c++
cmake -DBUILD_TESTING=ON ../src
```
<h4>Windows</h4>
The gtest library included in the repo needs to be built with 
forced shared libraries on Windows, so use the following:

```c++
cmake -DBUILD_TESTING=ON -Dgtest_forced_shared_crt=ON ../src
```
Be sure to include all other necessary CMake definitions as annotated above.
