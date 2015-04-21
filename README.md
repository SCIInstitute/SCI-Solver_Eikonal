SCI-Solver_Eikonal-3D
=====================

**Currently in pre-alpha stage, estimated stable release: Sep 2015**

SCI-Solver_Eikonal-3D is a C++/CUDA library written to solve the Eikonal equation on a 3D triangular mesh. It uses the fast iterative method (FIM) to solve efficiently, and uses GPU hardware.

The code was written by Zhisong Fu. The theory behind this code is published in the paper: "A Fast Iterative Method for Solving the Eikonal Equation on Tetrahedral Domains"

**AUTHORS:** Zhisong Fu(a,b), Robert M. Kirby(a,b), Ross T. Whitaker(a,b)

`  `a-School of Computing, University of Utah, Salt Lake City, UT, USA

`  `b-Scientific Computing and Imaging Institute, University of Utah, Salt Lake City, USA

**ABSTRACT:**
Generating numerical solutions to the eikonal equation and its many variations has a broad range of applications in both the natural and computational sciences. Efficient solvers on cutting-edge, parallel architectures require new algorithms that may not be theoretically optimal, but that are designed to allow asynchronous solution updates and have limited memory access patterns. This paper presents a parallel algorithm for solving the eikonal equation on fully unstructured tetrahedral meshes. The method is appropriate for the type of fine-grained parallelism found on modern massively-SIMD architectures such as graphics processors and takes into account the particular constraints and capabilities of these computing platforms. This work builds on previous work for solving these equations on triangle meshes; in this paper we adapt and extend previous 2D strategies to accommodate three-dimensional, unstructured, tetrahedralized domains. These new developments include a local update strategy with data compaction for tetrahedral meshes that provides solutions on both serial and parallel architectures, with a generalization to inhomogeneous, anisotropic speed functions. We also propose two new update schemes, specialized to mitigate the natural data increase observed when moving to three dimensions, and the data structures necessary for efficiently mapping data to parallel SIMD processors in a way that maintains computational density. Finally, we present descriptions of the implementations for a single CPU, as well as multicore CPUs with shared memory and SIMD architectures, with comparative results against state-of-the-art eikonal solvers.

Requirements
==============

 * Git, CMake (3.0+ recommended), and the standard system build environment tools.
 * You will need a CUDA Compatible Graphics card. See <a href="https://developer.nvidia.com/cuda-gpus">here</a> You will also need to be sure your card has CUDA compute capability of at least 2.0.
 * SCI-Solver_Eikonal-3D is compatible with the latest CUDA toolkit (7.0). Download <a href="https://developer.nvidia.com/cuda-downloads">here</a>.
 * This project has been tested on OpenSuse 12.3 (Dartmouth) on NVidia GeForce GTX 570 HD, Windows 7 on NVidia GeForce GTX 775M, and OSX 10.10 on NVidia GeForce GTX 775M. 
 * If you have a CUDA compatible card with the above operating systems, and are experiencing issues, please contact the repository owners.
 * Windows: You will need Microsoft Visual Studio 2013 build tools. This document describes the "NMake" process.
 * OSX: Please be sure to follow setup for CUDA <a href="http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-mac-os-x/#axzz3W4nXNNin">here</a>. There are several compatability requirements for different MAC machines, including using a different version of CUDA (ie. 5.5).

Building
==============

<h3>Unix / OSX</h3>
In a terminal:
```c++
mkdir SCI-SOLVER_Eikonal-3D/build
cd SCI-SOLVER_Eikonal-3D/build
cmake ../src
make
```

<h3>Windows</h3>
Open a Visual Studio (32 or 64 bit) Native Tools Command Prompt. 
Follow these commands:
```c++
mkdir C:\Path\To\SCI-Solver_Eikonal-3D\build
cd C:\Path\To\SCI-Solver_Eikonal-3D\build
cmake -G "NMake Makefiles" ..\src
nmake
```

**Note:** For all platforms, you may need to specify your CUDA toolkit location (especially if you have multiple CUDA versions installed):
```c++
cmake -DCUDA_TOOLKIT_ROOT_DIR="~/NVIDIA/CUDA-7.0" ../src
```
(Assuming this is the location).

**Note:** If you have compile errors such as <code>undefined reference: atomicAdd</code>, it is likely you need to set your compute capability manually. CMake outputs whether compute capability was determined automatically, or if you need to set it manually. The default (and known working) minimum compute capability is 2.0.

```c++
cmake -DCUDA_COMPUTE_CAPABILITY=20 ../src
make
```


Running Examples
==============

You will need to enable examples in your build to compile and run them

```c++
cmake -DBUILD_EIKONAL3D_SOLVER_EXAMPLES=ON ../src
make
```

You will find the example binaries built in the <code>build/examples</code> directory.

Run the examples in the build directory:

```c++
examples/Example1
examples/Example2
...
```

Follow the example source code in <code>src/examples</code> to learn how to use the library.

Using the Library
==============

A basic usage of the library links to the <code>libEikonal3D_CORE</code> library during build and 
includes the headers needed, which are usually no more than:

```c++
#include "setup_solver.h"
#include "cuda_resources.h"
```

Then a program would setup the FEM parameters using the 
<code>AMG_Config</code> object and call <code>setup_solver()</code> to generate
the answer matrices.

You will need to make sure your CMake/Makfile/Build setup knows where to point for the library and header files. See the examples and their CMakeLists.txt.

Testing
==============
The repo comes with a set of regression tests to see if recent changes break expected results. To build the tests, you will need to set <code>BUILD_TESTING</code> to "ON" in either <code>ccmake</code> or when calling CMake:

```c++
cmake -DBUILD_TESTING=ON ../src
```
<h4>Windows</h4>
The gtest library included in the repo needs to be built with forced shared libraries on Windows, so use the following:

```c++
cmake -DBUILD_TESTING=ON -Dgtest_forced_shared_crt=ON ../src
```
Be sure to include all other necessary CMake definitions as annotated above.
