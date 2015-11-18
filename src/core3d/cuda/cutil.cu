#include <cutil.h>
#include <cstdio>

void cudaSafeCall()
{
#ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err )
  {
    fprintf( stderr, "cudaSafeCall() failed at ");
    printFilenameAndLine();
    fprintf( stderr, " : %s\n", cudaGetErrorString( err ) );
    exit( -1 );
  }
#endif

  return;
}

void cudaSafeCall(cudaError err)
{
#ifdef CUDA_ERROR_CHECK
  if ( cudaSuccess != err )
  {
    fprintf( stderr, "cudaSafeCall() failed at ");
    printFilenameAndLine();
    fprintf( stderr, " : %s\n", cudaGetErrorString( err ) );
    exit( -1 );
  }
#endif

  return;
}
