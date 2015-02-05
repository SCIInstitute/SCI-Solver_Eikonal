SCI-Solver_Eikonal-3D
=====================

**Currently in pre-alpha stage, estimated stable release: March 2015**

SCI-Solver_Eikonal-3D is a C++/CUDA library written to solve the Eikonal equation on a 3D triangular mesh. It uses the fast iterative method (FIM) to solve efficiently, and uses GPU hardware.

The code was written by Zhisong Fu. The theory behind this code is published in the paper: "A Fast Iterative Method for Solving the Eikonal Equation on Tetrahedral Domains"

**AUTHORS:** Zhisong Fu(a,b), Robert M. Kirby(a,b), Ross T. Whitaker(a,b)

`  `a-School of Computing, University of Utah, Salt Lake City, UT, USA

`  `b-Scientific Computing and Imaging Institute, University of Utah, Salt Lake City, USA

**ABSTRACT:**
Generating numerical solutions to the eikonal equation and its many variations has a broad range of applications in both the natural and computational sciences. Efficient solvers on cutting-edge, parallel architectures require new algorithms that may not be theoretically optimal, but that are designed to allow asynchronous solution updates and have limited memory access patterns. This paper presents a parallel algorithm for solving the eikonal equation on fully unstructured tetrahedral meshes. The method is appropriate for the type of fine-grained parallelism found on modern massively-SIMD architectures such as graphics processors and takes into account the particular constraints and capabilities of these computing platforms. This work builds on previous work for solving these equations on triangle meshes; in this paper we adapt and extend previous 2D strategies to accommodate three-dimensional, unstructured, tetrahedralized domains. These new developments include a local update strategy with data compaction for tetrahedral meshes that provides solutions on both serial and parallel architectures, with a generalization to inhomogeneous, anisotropic speed functions. We also propose two new update schemes, specialized to mitigate the natural data increase observed when moving to three dimensions, and the data structures necessary for efficiently mapping data to parallel SIMD processors in a way that maintains computational density. Finally, we present descriptions of the implementations for a single CPU, as well as multicore CPUs with shared memory and SIMD architectures, with comparative results against state-of-the-art eikonal solvers.
