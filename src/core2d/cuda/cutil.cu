#include <cuda.h>
#include <cuda_runtime.h>
#include <cutil.h>
#include <thrust/scan.h>

template<typename IndexType, typename ValueType>
void __global__ neighbor_count_kernel(IndexType* tri0, IndexType* tri1, IndexType* tri2, IndexType ne, IndexType* nbcount)
{
  for(int eidx = threadIdx.x; eidx < ne; eidx += gridDim.x * blockDim.x)
  {
    IndexType i = tri0[eidx];
    IndexType j = tri1[eidx];
    IndexType k = tri2[eidx];

    atomicInc((unsigned *)nbcount + i, INT_MAX);
    atomicInc((unsigned *)nbcount + j, INT_MAX);
    atomicInc((unsigned *)nbcount + k, INT_MAX);

  }

}

template<typename IndexType, typename ValueType>
void __global__ compute_nb_indices_kernel(IndexType* rowoffsets, IndexType* ele_indices, IndexType *tri0, IndexType* tri1, IndexType* tri2, IndexType nv, IndexType* column_indices, size_t num_cols, size_t pitch)
{

  for(int nidx = threadIdx.x; nidx < nv; nidx += gridDim.x * blockDim.x)
  {
    for(int i = 0; i < num_cols; i++)
    {
      column_indices[pitch * i + nidx] = -1;
    }

    int nedges = 0;
    for(int j = rowoffsets[nidx]; j < rowoffsets[nidx + 1]; j++)
    {
      IndexType jj = ele_indices[j];
      IndexType node0 = tri0[jj];
      IndexType node1 = tri1[jj];
      IndexType node2 = tri2[jj];
      if(node0 != nidx)
      {
        column_indices[pitch * nedges + nidx] = node0;
        nedges++;
      }

    }
  }

}

template<typename IndexType, typename ValueType>
void __global__ compute_ele_indices_kernel(IndexType* tri0, IndexType* tri1, IndexType* tri2, IndexType ne, IndexType* rowoffsets, IndexType* ele_indices)
{
  for(int eidx = threadIdx.x; eidx < ne; eidx += gridDim.x * blockDim.x)
  {
    IndexType i = tri0[eidx];
    IndexType j = tri1[eidx];
    IndexType k = tri2[eidx];
    IndexType starti = rowoffsets[i];
    IndexType startj = rowoffsets[j];
    IndexType startk = rowoffsets[k];
    IndexType endi = rowoffsets[i + 1];
    IndexType endj = rowoffsets[j + 1];
    IndexType endk = rowoffsets[k + 1];
    for(int n = starti; n < endi; n++)
    {
      atomicCAS(ele_indices + n, -1, eidx);
      break;
    }

    for(int n = startj; n < endj; n++)
    {
      atomicCAS(ele_indices + n, -1, eidx);
      break;
    }

    for(int n = startk; n < endk; n++)
    {
      atomicCAS(ele_indices + n, -1, eidx);
    }


  }

}

template<typename IndexType, typename ValueType>
__global__ void convert_kernel(IndexType* rowoff1, IndexType* colidx1, ValueType* values1, IndexType* rowidx2, IndexType* colidx2, ValueType* values2, int num_rows)
{
  for(int ridx = blockIdx.x * blockDim.x + threadIdx.x; ridx < num_rows; ridx++)
  {
    IndexType start1 = rowoff1[ridx];
    IndexType end1 = rowoff1[ridx + 1];
    IndexType start2 = start1 * 2 - ridx;

    rowidx2[start2] = ridx;
    colidx2[start2] = ridx;
    values2[start2] = values1[start1];
    for(int i = start1 + 1; i < end1; i++)
    {
      ValueType v = values1[i];
      IndexType col = colidx1[i];
      IndexType loc = start2 + 1 + 2 * (i - start1 - 1);
      rowidx2[loc] = ridx;
      colidx2[loc] = col;
      values2[loc] = v;
      rowidx2[loc + 1] = col;
      colidx2[loc + 1] = ridx;
      values2[loc + 1] = v;
    }
  }

}
