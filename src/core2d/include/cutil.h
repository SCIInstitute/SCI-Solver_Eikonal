#ifndef __CUTIL_H__
#define __CUTIL_H__

#include <TriMesh.h>



/**********************************************************
 * Checks for a cuda error and if one exists prints it,
 * the stack trace, and exits
 *********************************************************/
#define CUDA_ERROR_CHECK

#define cudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define cudaCheckErrors()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
  if ( cudaSuccess != err )
  {
    fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
        file, line, cudaGetErrorString( err ) );
    exit( -1 );
  }
#endif

  return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err )
  {
    fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
        file, line, cudaGetErrorString( err ) );
    exit( -1 );
  }

  // More careful checking. However, this will affect performance.
  // Comment away if needed.
  err = cudaDeviceSynchronize();
  if( cudaSuccess != err )
  {
    fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
        file, line, cudaGetErrorString( err ) );
    exit( -1 );
  }
#endif

  return;
}

template <class Matrix, class Vector>
void computeResidual(const Matrix& A, const Vector& x, const Vector& b, Vector& r);

template<typename IndexType, typename ValueType>
__global__ void find_diag_kernel(const IndexType num_rows, const IndexType num_cols, const IndexType num_cols_per_row, const IndexType pitch,
                                 const IndexType * Aj,
                                 const ValueType* Ax,
                                 ValueType* diag)
{
  const IndexType thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  const IndexType grid_size = gridDim.x * blockDim.x;

  for (IndexType row = thread_id; row < num_rows; row += grid_size)
  {
    IndexType offset = row;

    for (IndexType n = 0; n < num_cols_per_row; n++)
    {
      const IndexType col = Aj[offset];

      if (col == row)
      {
        const ValueType A_ij = Ax[offset];
        diag[row] = A_ij;
      }

      offset += pitch;
    }
  }
}

/**************************************************
 * structs for converting between signed and unsigned values without 
 * type casting.
 * ************************************************/

/*****************************
 * Generic converter for unsigned types.
 * This becomes a no op
 *****************************/
template <class GlobalOrdinal>
struct intuint
{

  union
  {
    GlobalOrdinal ival;
    GlobalOrdinal uval;
  };
};

/***************************
 * char converter
 **************************/
template <>
struct intuint<char>
{

  union
  {
    char ival;
    unsigned char uval;
  };
};

/***************************
 * Short converter
 **************************/
template <>
struct intuint<short>
{

  union
  {
    short ival;
    unsigned short uval;
  };
};

/***************************
 * Integer converter
 **************************/
template <>
struct intuint<int>
{

  union
  {
    int ival;
    unsigned int uval;
  };
};

/***************************
 * long converter
 **************************/
template <>
struct intuint<long>
{

  union
  {
    long ival;
    unsigned long uval;
  };
};

struct metisinput
{
  int nn;
  int* xadj;
  int* adjncy;
};

struct cudaCSRGraph
{
  int nn;
  int* xadj;
  int* adjncy;
};

#endif
