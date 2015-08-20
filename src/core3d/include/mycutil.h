#ifndef MYCUTIL_H
#define MYCUTIL_H

//#include <myerror.h>
/**********************************************************
 * Checks for a cuda error and if one exists prints it,
 * the stack trace, and exits
 *********************************************************/
#define FatalError(x)
#define cudaCheckError() {                              \
  cudaError_t e=cudaGetLastError();                                 \
  char error_str[100];                                              \
  if(e!=cudaSuccess) {                                              \
    sprintf(error_str,"Cuda failure: '%s'",cudaGetErrorString(e));  \
    FatalError(error_str);                                          \
  }                                                                 \
}                                                                   \

#define cudaSafeCall(x) {(x); cudaCheckError()}


template <typename T>
__inline__ __device__ __host__ T DOT_PRODUCT(const T* a, const T* b)
{
  T c = 0;
  for (int i = 0; i < 3; i++)
    c += a[i] * b[i];
  return c;
}

template <typename T>
__inline__ __device__ __host__ void CROSS_PRODUCT(const T* v1, const T* v2, T* v3)
{
  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

template <typename T>
__inline__ __device__ __host__ T LENGTH(const T* v)
{
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}


#endif
