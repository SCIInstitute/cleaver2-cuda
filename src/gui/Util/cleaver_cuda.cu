/**
 * @file cleaver_cuda.cu
 * This file contains the GPU calls for Cleaver.
 * @version April 21, 2014
 * @author: Brig Bagley
 */
#include <cstdio>
#include <iostream>
#include <stdint.h>
#include "cleaver_cuda.hh"

namespace CleaverCUDA {
/**
 * The error checker for CUDA.
 * @param val The return value from a CUDA Call.
 */
void CudaCheckReturn(cudaError_t val) {
  if (val != cudaSuccess) {
    fprintf(stderr, "Error %s at line %d in file %s\n",
            cudaGetErrorString(val), __LINE__, __FILE__);
    exit(1);
  }
}

__device__ __host__
float DataTransformCUDA(float *data_,
                        float i, float j, float k,
                        size_t mm, size_t m, float* scales,
                        float* scale,
                        size_t ww_, size_t hh_, size_t dd_) {

  float x = i;
  float y = j;
  float z = k;

  size_t whd = ww_*hh_*dd_;
  size_t wh = ww_*hh_;
  size_t w = ww_;
  size_t h = hh_;
  size_t d = dd_;

  x *= scale[0];
  y *= scale[1];
  z *= scale[2];
  if (mm < m - 1) {
    x *= scales[mm*3+0];
    y *= scales[mm*3+1];
    z *= scales[mm*3+2];
  }

  x -= 0.5f;
  y -= 0.5f;
  z -= 0.5f;

  bool inside = (x >= -.5) && (x < ww_+.5) &&
      (y >= -.5) && (y < hh_+.5) &&
      (z >= -.5) && (z < dd_+.5);

  if (mm < m - 1) {
    if(inside) {
      float t = fmodf(x,1.0f);
      float u = fmodf(y,1.0f);
      float v = fmodf(z,1.0f);

      int i0 = (int)(floorf(x));   int i1 = i0+1;
      int j0 = (int)(floorf(y));   int j1 = j0+1;
      int k0 = (int)(floorf(z));   int k1 = k0+1;
      int zero = 0;

      i0 = min(max(zero,i0), (int)(w)-1);
      j0 = min(max(zero,j0),(int)(h)-1);
      k0 = min(max(zero,k0),(int)(d)-1);

      i1 = min(max(zero,i1),(int)(w)-1);
      j1 = min(max(zero,j1),(int)(h)-1);
      k1 = min(max(zero,k1),(int)(d)-1);

      float C000 = data_[i0 + j0*w + k0*wh + mm*whd];
      float C001 = data_[i0 + j0*w + k1*wh + mm*whd];
      float C010 = data_[i0 + j1*w + k0*wh + mm*whd];
      float C011 = data_[i0 + j1*w + k1*wh + mm*whd];
      float C100 = data_[i1 + j0*w + k0*wh + mm*whd];
      float C101 = data_[i1 + j0*w + k1*wh + mm*whd];
      float C110 = data_[i1 + j1*w + k0*wh + mm*whd];
      float C111 = data_[i1 + j1*w + k1*wh + mm*whd];

      return float((1-t)*(1-u)*(1-v)*C000 + (1-t)*(1-u)*(v)*C001 +
                   (1-t)*  (u)*(1-v)*C010 + (1-t)*  (u)*(v)*C011 +
                   (t)*(1-u)*(1-v)*C100 +   (t)*(1-u)*(v)*C101 +
                   (t)*  (u)*(1-v)*C110 +   (t)*  (u)*(v)*C111);
    } else
      return -1000.;
  } else {
    if(inside)
      return -1000.;
    else
      return 1000.;
  }

}

/**
 * Runs the loops of chunks to find max materials.
 * @param input_device_memory Pointer to the material data.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @param m The input data number of materials.
 * @param output_device_memory The pointer to the chunk of results.
 * @param i The x location of the start of the chunk.
 * @param j The y location of the start of the chunk.
 * @param k The z location of the start of the chunk.
 * @param endi The x location of the end of the chunk.
 * @param endj The y location of the end of the chunk.
 * @param endk The z location of the end of the chunk.
 * @param device_scales The pointer to the scales to use.
 * @param device_scale The pointer to the scale to use.
 */
__global__
void FindMaxesCUDA(void *input_device_memory,
                   size_t w, size_t h, size_t d, size_t m,
                   void *output_device_memory,
                   size_t i, size_t j, size_t k,
                   size_t endi, size_t endj, size_t endk,
                   void *device_scales, void * device_scale) {
  //cast to proper data
  float* data_ = (float*)input_device_memory;
  float* scales = (float*)device_scales;
  float* scale = (float*)device_scale;
  char* labels = (char*)output_device_memory;

  size_t grid_max = gridDim.x * gridDim.y * gridDim.z;

  size_t idx =
      blockIdx.x       +
      blockIdx.y * gridDim.x  +
      blockIdx.z * gridDim.x * gridDim.y  +
      threadIdx.x  * grid_max  +
      threadIdx.y  * grid_max * blockDim.x +
      threadIdx.z  * grid_max * blockDim.x * blockDim.y;

  size_t x = idx % kChunkSize;
  idx -= x; idx /= kChunkSize;
  size_t y = idx % kChunkSize;
  idx -= y; idx /= kChunkSize;
  size_t z = idx % kChunkSize;

  //  for (size_t x = 0; x < (kBlockSize); x++ )
  //    for (size_t y = 0; y < (kBlockSize); y++ )
  //      for (size_t z = 0; z < (kBlockSize); z++ )
  if ((x < endi - i) && (y < endj - j) && (z < endk - k))
    labels[x +
           y * kChunkSize +
           z * kChunkSize * kChunkSize] =
               GetLabelCUDA(data_,scales, scale,i+x,j+y,k+z,m,w,h,d);
}

__device__ __host__
float GetCellCenterValueCUDA(
    float* data_, float* scales, float* scale,
    size_t i, size_t j, size_t k,
    size_t m, size_t num_mats,
    bool find_max, char* max_mat,
    size_t w, size_t h, size_t d) {
  if (!find_max) {
    *max_mat = m;
    return DataTransformCUDA(
        data_,
        static_cast<float>(i)+0.5f,
        static_cast<float>(j)+0.5f,
        static_cast<float>(k)+0.5f,
        m,num_mats,scales,scale,w,h,d);
  }
  char mat = 0;
  float tmp, val = DataTransformCUDA(
      data_,
      static_cast<float>(i)+0.5f,
      static_cast<float>(j)+0.5f,
      static_cast<float>(k)+0.5f,
      0,num_mats,scales,scale,w,h,d);
  for(size_t a = 1; a < num_mats; a++) {
    if ((tmp = DataTransformCUDA(
        data_,
        static_cast<float>(i)+0.5f,
        static_cast<float>(j)+0.5f,
        static_cast<float>(k)+0.5f,
        a,num_mats,scales,scale,w,h,d)) > val) {
      val = tmp;
      mat = a;
    }
  }
  *max_mat = mat;
  return val;
}

__device__ __host__
void SetArrayCUDA(size_t *cell, size_t i, size_t j, size_t k) {
  cell[0] = i;
  cell[1] = j;
  cell[2] = k;
}

__device__ __host__
void GetAdjacentCellFromEdgeCUDA(
    CleaverCUDA::edge_index num, bool first,
    size_t i, size_t j, size_t k, size_t* cell) {
  switch (num) {
    //Diagonal edges
    case CleaverCUDA::DULF:
      SetArrayCUDA(cell,i,j+1,k); return;
    case CleaverCUDA::DULB:
      SetArrayCUDA(cell,i,j+1,k+1); return;
    case CleaverCUDA::DURF:
      SetArrayCUDA(cell,i+1,j+1,k); return;
    case CleaverCUDA::DURB:
      SetArrayCUDA(cell,i+1,j+1,k+1); return;
    case CleaverCUDA::DLLF:
      SetArrayCUDA(cell,i,j,k); return;
    case CleaverCUDA::DLLB:
      SetArrayCUDA(cell,i,j,k+1); return;
    case CleaverCUDA::DLRF:
      SetArrayCUDA(cell,i+1,j,k); return;
    case CleaverCUDA::DLRB:
      SetArrayCUDA(cell,i+1,j,k+1); return;
      //Dual Edges
    case CleaverCUDA::CL:
      SetArrayCUDA(cell,i-1,j,k); return;
    case CleaverCUDA::CR:
      SetArrayCUDA(cell,i+1,j,k); return;
    case CleaverCUDA::CU:
      SetArrayCUDA(cell,i,j+1,k); return;
    case CleaverCUDA::CD:
      SetArrayCUDA(cell,i,j-1,k); return;
    case CleaverCUDA::CF:
      SetArrayCUDA(cell,i,j,k-1); return;
    case CleaverCUDA::CB:
      SetArrayCUDA(cell,i,j,k+1); return;
      //Axis edges (top)
    case CleaverCUDA::UL:
      if (first) SetArrayCUDA(cell,i,j+1,k);
      else
        SetArrayCUDA(cell,i,j+1,k+1);
      return;
    case CleaverCUDA::UR:
      if (first) SetArrayCUDA(cell, i+1,j+1,k);
      else
        SetArrayCUDA(cell,i+1,j+1,k+1);
      return;
    case CleaverCUDA::UF:
      if (first) SetArrayCUDA(cell,i,j+1,k);
      else
        SetArrayCUDA(cell,i+1,j+1,k);
      return;
    case CleaverCUDA::UB:
      if (first) SetArrayCUDA(cell,i,j+1,k+1);
      else
        SetArrayCUDA(cell,i+1,j+1,k+1);
      return;
      //axis edges (bottom)
    case CleaverCUDA::LL:
      if (first) SetArrayCUDA(cell,i,j,k);
      else
        SetArrayCUDA(cell,i,j,k+1);
      return;
    case CleaverCUDA::LR:
      if (first) SetArrayCUDA(cell,i+1,j,k);
      else
        SetArrayCUDA(cell,i+1,j,k+1);
      return;
    case CleaverCUDA::LF:
      if (first) SetArrayCUDA(cell,i,j,k);
      else
        SetArrayCUDA(cell,i+1,j,k);
      return;
    case CleaverCUDA::LB:
      if (first) SetArrayCUDA(cell,i,j,k+1);
      else
        SetArrayCUDA(cell,i+1,j,k+1);
      return;
      //axis edges (columns)
    case CleaverCUDA::FL:
      if (first) SetArrayCUDA(cell,i,j,k);
      else
        SetArrayCUDA(cell,i,j+1,k);
      return;
    case CleaverCUDA::FR:
      if (first) SetArrayCUDA(cell,i+1,j,k);
      else
        SetArrayCUDA(cell,i+1,j+1,k);
      return;
    case CleaverCUDA::BL:
      if (first) SetArrayCUDA(cell,i,j,k+1);
      else
        SetArrayCUDA(cell,i,j+1,k+1);
      return;
    case CleaverCUDA::BR:
      if (first) SetArrayCUDA(cell,i+1,j,k+1);
      else
        SetArrayCUDA(cell,i+1,j+1,k+1);
      return;
  }
}

__device__ __host__
char GetLabelCUDA(float* data_, float* scales, float* scale,
              size_t i, size_t j, size_t k, char num_mats,
              size_t w, size_t h, size_t d) {
  float max = DataTransformCUDA(data_, static_cast<float>(i),
                                static_cast<float>(j),
                                static_cast<float>(k),
                                0,num_mats,scales,scale,w,h,d);
  char max_mat = 0;
  for(char t = 1; t < num_mats; t++) {
    float tmp;
    if ((tmp = DataTransformCUDA(data_, static_cast<float>(i),
                                 static_cast<float>(j),
                                 static_cast<float>(k),
                                 t,num_mats,scales,scale,w,h,d)) > max) {
      max = tmp;
      max_mat = t;
    }
  }
  return max_mat;
}

__device__ __host__
float GetEdgeMatAndValueAtEndpointCUDA(
    float* data_, float* scales, float* scale,
    CleaverCUDA::edge_index num, bool first, bool find_max,
    size_t m, size_t num_mats, char* ret_mat,
    size_t i, size_t j, size_t k,
    size_t w, size_t h, size_t d) {
  //find the respective adjacent cell for the vertex we want
  size_t arr[3];
  CleaverCUDA::GetAdjacentCellFromEdgeCUDA(num,first,i,j,k,arr);
  //the center value of this cell. (V1)
  char mat;
  float res = CleaverCUDA::GetCellCenterValueCUDA(
      data_,scales,scale,i,j,k,
      m,num_mats,find_max,&mat,w,h,d);
  if ((num < 8 && !first) || (num >= 14)) {
    // diagonal edges & second vertex, or axis edges.
    mat = find_max?CleaverCUDA::GetLabelCUDA(
        data_,scales,scale,arr[0],arr[1],arr[2],
        num_mats,w,h,d):m;
    res = CleaverCUDA::DataTransformCUDA(
        data_,arr[0],arr[1],arr[2],mat,num_mats,
        scales,scale,w,h,d);
  } else if ((8 <= num && num < 14)) { // dual edges
    bool neg =
        (num == CleaverCUDA::CL) || // is always on left, bottom, or front.
        (num == CleaverCUDA::CD) ||
        (num == CleaverCUDA::CF);
    if (neg) first = !first;       // flip first if we are on the negative edge
    if (!first)                    // second vertex
      res = CleaverCUDA::GetCellCenterValueCUDA(
          data_,scales,scale,arr[0],arr[1],arr[2],
          m,num_mats,find_max,&mat,w,h,d);
  }
  *ret_mat = mat;
  return res;
}

__device__ __host__
void GetEdgeVerticesCUDA(
    CleaverCUDA::edge_index num,
    size_t i, size_t j, size_t k,
    float* scale, float verts[2][3]) {
  float x = static_cast<float>(i) * scale[0];
  float y = static_cast<float>(j) * scale[1];
  float z = static_cast<float>(k) * scale[2];
  float v1[3] = {x,y,z};
  float v2[3] = {x+scale[0],y,z};
  float v3[3] = {x,y+scale[1],z};
  float v4[3] = {x+scale[0],y+scale[1],z};
  float v5[3] = {x,y,z+scale[2]};
  float v6[3] = {x+scale[0],y,z+scale[2]};
  float v7[3] = {x,y+scale[1],z+scale[2]};
  float v8[3] = {x+scale[0],y+scale[1],z+scale[2]};
  float v9[3] = {x+.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]};
  float v10[3] = {x-.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]};
  float v11[3] = {x+1.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]};
  float v12[3] = {x+.5f * scale[0],
      y-.5f * scale[1],z+.5f * scale[2]};
  float v13[3] = {x+.5f * scale[0],
      y+1.5f * scale[1],z+.5f * scale[2]};
  float v14[3] = {x+.5f * scale[0],
      y+.5f * scale[1],z-.5f * scale[2]};
  float v15[3] = {x+.5f * scale[0],
      y+.5f * scale[1],z+1.5f * scale[2]};

  float* ans[2] = {NULL, NULL};
  switch (num) {
    //diagonal edges upper
    case CleaverCUDA::DULF: ans[0] = v9; ans[1] = v3; break;
    case CleaverCUDA::DULB: ans[0] = v9; ans[1] = v7; break;
    case CleaverCUDA::DURF: ans[0] = v9; ans[1] = v4; break;
    case CleaverCUDA::DURB: ans[0] = v9; ans[1] = v8; break;
    //diagonal edges lower
    case CleaverCUDA::DLLF: ans[0] = v9; ans[1] = v1; break;
    case CleaverCUDA::DLLB: ans[0] = v9; ans[1] = v5; break;
    case CleaverCUDA::DLRF: ans[0] = v9; ans[1] = v2; break;
    case CleaverCUDA::DLRB: ans[0] = v9; ans[1] = v6; break;
    //dual edges
    case CleaverCUDA::CL: ans[0] = v10; ans[1] = v9; break;
    case CleaverCUDA::CR: ans[0] = v9; ans[1] = v11; break;
    case CleaverCUDA::CU: ans[0] = v9; ans[1] = v13; break;
    case CleaverCUDA::CD: ans[0] = v12; ans[1] = v9; break;
    case CleaverCUDA::CF: ans[0] = v14; ans[1] = v9; break;
    case CleaverCUDA::CB: ans[0] = v9; ans[1] = v15; break;
    //top face edges
    case CleaverCUDA::UL: ans[0] = v3; ans[1] = v7; break;
    case CleaverCUDA::UR: ans[0] = v4; ans[1] = v8; break;
    case CleaverCUDA::UF: ans[0] = v3; ans[1] = v4; break;
    case CleaverCUDA::UB: ans[0] = v7; ans[1] = v8; break;
    //bottom face edges
    case CleaverCUDA::LL: ans[0] = v1; ans[1] = v5; break;
    case CleaverCUDA::LR: ans[0] = v2; ans[1] = v6; break;
    case CleaverCUDA::LF: ans[0] = v1; ans[1] = v2; break;
    case CleaverCUDA::LB: ans[0] = v5; ans[1] = v6; break;
    //column edges
    case CleaverCUDA::FL: ans[0] = v1; ans[1] = v3; break;
    case CleaverCUDA::FR: ans[0] = v2; ans[1] = v4; break;
    case CleaverCUDA::BL: ans[0] = v5; ans[1] = v7; break;
    case CleaverCUDA::BR: ans[0] = v6; ans[1] = v8; break;
  }
  for (size_t t = 0; t < 2; t++)
    for (size_t tt = 0; tt < 3; tt++)
      verts[t][tt] = ans[t][tt];
}

__device__ __host__
void FindEdgeCutCUDA(
    float* data_,
    float* scales,
    float* scale,
    size_t w, size_t h, size_t d, size_t m,
    size_t i, size_t j, size_t k,
    Edge* edge, CleaverCUDA::edge_index num) {
  edge->isCut_eval |= CleaverCUDA::kIsEvaluated;
  //get strongest material at each end of the edge
  unsigned char matA, matB, dummy;
  float v1 = CleaverCUDA::GetEdgeMatAndValueAtEndpointCUDA(
      data_,scales,scale,num,true,true,
      0,m,(char*)&matA,i,j,k,w,h,d);
  float v2 = CleaverCUDA::GetEdgeMatAndValueAtEndpointCUDA(
      data_,scales,scale,num,false,true,
      0,m,(char*)&matB,i,j,k,w,h,d);
  //if they are the same, nothing to be done: no cut
  if (matA == matB) return;
  //if they are different, interpolate transition point.
  float a1 = v1;
  float b2 = v2;
  float a2 = CleaverCUDA::GetEdgeMatAndValueAtEndpointCUDA(
      data_,scales,scale,num,false,false,
      matA,m,(char*)&dummy,i,j,k,w,h,d);
  float b1 = CleaverCUDA::GetEdgeMatAndValueAtEndpointCUDA(
      data_,scales,scale,num,true,false,
      matB,m,(char*)&dummy,i,j,k,w,h,d);
  float top = (a1 - b1);
  float bot = (b2 - a2 + a1 - b1);
  //degenerate cases
  if (bot == 0.) return;
  edge->isCut_eval |= CleaverCUDA::kIsCut;
  matA = matA % kMaxMaterials;
  matB = matB % kMaxMaterials;
  unsigned char minMat = min(matA,matB);
  unsigned char maxMat = max(matA,matB);
  unsigned char adder = 0;
  switch(minMat) {
    case 1: adder = 4; break;
    case 2: adder = 7; break;
    case 3: adder = 9; break;
    case 4: adder = 10; break;
    default: adder = 0; break;
  }
  unsigned char mm = (maxMat - 1 + adder)%16;
  edge->isCut_eval |= (mm << CleaverCUDA::kMaterial);
  float t = min(max(top/bot,0.f),1.f);
  float edge_verts[2][3];
  CleaverCUDA::GetEdgeVerticesCUDA(num,i,j,k,scale,edge_verts);
  for (size_t x = 0; x < 3; x++)
    edge->cut_loc[x] = (1. - t) * edge_verts[0][x] + t * edge_verts[1][x];
}

/**
 * This function determines all of the inner edge cuts.
 * @param input_device_memory Pointer to the material data.
 * @param device_scales The pointer to the scales to use.
 * @param device_scale The pointer to the scale to use.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @param m The input data number of materials.
 * @param i The x location of the start of the chunk.
 * @param j The y location of the start of the chunk.
 * @param k The z location of the start of the chunk.
 * @param endi The x location of the end of the chunk.
 * @param endj The y location of the end of the chunk.
 * @param endk The z location of the end of the chunk.
 * @param output_device_memory The pointer to the chunk of results.
 * @param which Which set of edges we're working on (inner, dual, or axis).
 */
__global__
void FindEdgeCutsCUDA(
    void* input_device_memory,
    void* device_scales,
    void* device_scale,
    size_t w, size_t h, size_t d, size_t m,
    size_t i, size_t j, size_t k,
    size_t endi, size_t endj, size_t endk,
    void* output_device_memory,
    char which) {
  //cast to proper data
  float* data_ = (float*)input_device_memory;
  float* scales = (float*)device_scales;
  float* scale = (float*)device_scale;
  Edge* edges = (Edge*)output_device_memory;
  size_t grid_max = gridDim.x * gridDim.y * gridDim.z;

  size_t idx =
      blockIdx.x       +
      blockIdx.y * gridDim.x  +
      blockIdx.z * gridDim.x * gridDim.y  +
      threadIdx.x  * grid_max  +
      threadIdx.y  * grid_max * blockDim.x +
      threadIdx.z  * grid_max * blockDim.x * blockDim.y;

  size_t num_edges = (which=='i')?8:3;

  size_t e = idx % num_edges;
  idx -= e; idx /= num_edges;
  size_t x = idx % kChunkSize;
  idx -= x; idx /= kChunkSize;
  size_t y = idx % kChunkSize;
  idx -= y; idx /= kChunkSize;
  size_t z = idx % kChunkSize;

  CleaverCUDA::edge_index edge_nums[8] =
  {DULF, DULB, DURF, DURB, DLLF, DLLB, DLRF, DLRB};
  if (which=='d') {
    edge_nums[0] = CL;
    edge_nums[1] = CD;
    edge_nums[2] = CF;
  } else if (which == 'a') {
    edge_nums[0] = LF;
    edge_nums[1] = FL;
    edge_nums[2] = LL;
  }
  //    for (size_t x = 0; x < kBlockSize; x++ )
  //      for (size_t y = 0; y < kBlockSize; y++ )
  //        for (size_t z = 0; z < kBlockSize; z++ ) {
//  for (size_t e = 0; e < num_edges; e++) {
    CleaverCUDA::edge_index edge_num = edge_nums[e];
    Edge *edge = &edges[(
        x +
        y * kChunkSize +
        z * kChunkSize * kChunkSize) * num_edges + e];
    //first clear the data
    edge->isCut_eval = 0;
    edge->cut_loc[0] = edge->cut_loc[1] = edge->cut_loc[2] = 0.0;
    if ((x < endi - i) && (y < endj - j) &&
        (z < endk - k) && (e < num_edges)) {
      FindEdgeCutCUDA(data_,scales,scale,w,h,d,m,
                      i+x,j+y,k+z,edge,edge_num);
    }
//  }
  //        }
}

__host__
void CallCUDAMaxes(float *all_data,
                   float * scales, float * scale,
                   size_t w, size_t h, size_t d, size_t m,
                   char * labels,
                   size_t wl, size_t hl, size_t dl,
                   void* device_pointers[3]) {
  //The chunk size will always be static
  size_t chunk_size = kChunkSize;
  //declare and allocate device memory for NRRD input data.
  size_t input_memory_size = w * h * d * (m - 1);
  void *input_device_memory = NULL;
  CudaCheckReturn(cudaMalloc((void**)&input_device_memory,
                             sizeof(float) *
                             input_memory_size));
  CudaCheckReturn(cudaMemcpy(input_device_memory, all_data,
                             sizeof(float) *
                             input_memory_size,
                             cudaMemcpyHostToDevice));

  //declare and allocate memory for the output of the device.
  size_t output_memory_size = chunk_size * chunk_size * chunk_size;
  void *output_device_memory = NULL;
  CudaCheckReturn(cudaMalloc((void**)&output_device_memory,
                             sizeof(char) *
                             output_memory_size));
  // allocate, set up , and copy scales.
  void *device_scales = NULL;
  CudaCheckReturn(cudaMalloc((void**)&device_scales,
                             sizeof(float) * (m - 1) * 3));
  CudaCheckReturn(cudaMemcpy(device_scales, scales,
                             sizeof(float) * (m - 1) * 3,
                             cudaMemcpyHostToDevice));
  void *device_scale = NULL;
  CudaCheckReturn(cudaMalloc((void**)&device_scale,
                             sizeof(float) * 3));
  CudaCheckReturn(cudaMemcpy(device_scale, scale,
                             sizeof(float) * 3,
                             cudaMemcpyHostToDevice));
  //set up blocks and grid
  dim3 dimBlock(kThreadSize, kThreadSize, kThreadSize);
  dim3 dimGrid(kBlockSize, kBlockSize, kBlockSize);
  //for each block
  for (size_t i = 0;; i+=chunk_size) {
    size_t endi = i + chunk_size;
    if (i >= wl) break;
    if (endi > wl)
      endi = wl;
    if (endi - 1 == i) break;
    for (size_t j = 0;; j+=chunk_size) {
      size_t endj = j + chunk_size;
      if (j >= hl) break;
      if (endj > hl)
        endj = hl;
      if (endj - 1 == j) break;
      for (size_t k = 0;; k+=chunk_size) {
        size_t endk = k + chunk_size;
        if (k >= dl) break;
        if (endk > dl)
          endk = dl;
        if (endk - 1 == k) break;
        //call the kernel
        FindMaxesCUDA<<<dimGrid,dimBlock>>>(
            input_device_memory, w, h, d, m,
            output_device_memory, i, j, k, endi, endj, endk,
            device_scales, device_scale);
        CudaCheckReturn(cudaThreadSynchronize());
        //copy the block results back
        char lbls[output_memory_size];
        //        FindMaxesCUDA(
        //            all_data, w, h, d, m,
        //            lbls, i, j, k, endi, endj, endk,
        //            scales, scale);
        CudaCheckReturn(cudaMemcpy(lbls,output_device_memory,
                                   sizeof(char) *
                                   output_memory_size ,
                                   cudaMemcpyDeviceToHost));
        for (size_t ii = 0; ii < endi - i; ii++)
          for (size_t jj = 0; jj < endj - j; jj++)
            for (size_t kk = 0; kk < endk - k; kk++) {
              labels[(i + ii) + (j + jj)*wl + (k + kk)*wl*hl] =
                  lbls[ii + jj * chunk_size + kk * chunk_size * chunk_size];
            }
      }
    }
  }
  // Free no longer used GPU memory
  CudaCheckReturn(cudaFree(output_device_memory));
  device_pointers[0] = input_device_memory;
  device_pointers[1] = device_scales;
  device_pointers[2] = device_scale;
}

__host__
size_t CallCUDACuts(void* data,
                    void* scales, void* scale,
                    size_t w, size_t h, size_t d, size_t m,
                    size_t wl, size_t hl, size_t dl,
                    CleaverCUDA::Edge* inner_edges,
                    CleaverCUDA::Edge* dual_edges,
                    CleaverCUDA::Edge* axis_edges,
                    bool* cut_cells) {
  //The chunk size will always be static
  size_t chunk = kChunkSize;
  //declare and allocate output memory
  size_t output_memory_size = (chunk + 1) * (chunk + 1) * (chunk + 1) * 8;
  void *output_device_memory = NULL;
  CudaCheckReturn(cudaMalloc((void**)&output_device_memory,
                             sizeof(Edge) *
                             output_memory_size));
  //the number of cuts found
  size_t count = 0, max_cell = wl*hl*dl;
  //set up blocks and grid
  dim3 dimBlock(kThreadSize, kThreadSize,kThreadSize);
  dim3 dimGrid8(kBlockSize*2, kBlockSize*2,kBlockSize*2);
  dim3 dimGrid3(kBlockSize*3, kBlockSize,kBlockSize);
  CleaverCUDA::Edge *inner_edges_output =
      new CleaverCUDA::Edge[output_memory_size];
  CleaverCUDA::Edge *dual_edges_output =
      new CleaverCUDA::Edge[output_memory_size];
  CleaverCUDA::Edge *axis_edges_output =
      new CleaverCUDA::Edge[output_memory_size];
  for (size_t t = 0; t < output_memory_size; t++) {
    inner_edges_output[t].isCut_eval = 0;
    dual_edges_output[t].isCut_eval = 0;
    axis_edges_output[t].isCut_eval = 0;
    inner_edges_output[t].cut_loc[0] =
        inner_edges_output[t].cut_loc[1] =
            inner_edges_output[t].cut_loc[2] = 0.0;
    dual_edges_output[t].cut_loc[0] =
        dual_edges_output[t].cut_loc[1] =
            dual_edges_output[t].cut_loc[2] = 0.0;
    axis_edges_output[t].cut_loc[0] =
        axis_edges_output[t].cut_loc[1] =
            axis_edges_output[t].cut_loc[2] = 0.0;
  }
  //for each block
  for (size_t i = 0;; i+=chunk) {
    size_t endi = i + chunk;
    if (i >= wl) break;
    if (endi > wl)
      endi = wl;
    if (endi - 1 == i) break;
    for (size_t j = 0;; j+=chunk) {
      size_t endj = j + chunk;
      if (j >= hl) break;
      if (endj > hl)
        endj = hl;
      if (endj - 1 == j) break;
      for (size_t k = 0;; k+=chunk) {
        size_t endk = k + chunk;
        if (k >= dl) break;
        if (endk > dl)
          endk = dl;
        if (endk - 1 == k) break;
        //call the inner edge kernel
        FindEdgeCutsCUDA<<<dimGrid8,dimBlock>>>(
            data,
            scales,
            scale,
            w, h, d, m,
            i, j, k, endi, endj, endk,
            output_device_memory,'i');
        //                FindEdgeCutsCUDA(
        //                    input_device_memoryV,
        //                    device_scalesV,
        //                    device_scaleV,
        //                    w, h, d, m,
        //                    i, j, k, endi, endj, endk,
        //                    inner_edges_output,'i');
        CudaCheckReturn(cudaThreadSynchronize());
        //copy the block results back
        CudaCheckReturn(cudaMemcpy(inner_edges_output,
                                   output_device_memory,
                                   sizeof(Edge) *
                                   output_memory_size ,
                                   cudaMemcpyDeviceToHost));
        //call the dual edge kernel
        FindEdgeCutsCUDA<<<dimGrid3,dimBlock>>>(
            data,
            scales,
            scale,
            w, h, d, m,
            i, j, k, endi, endj, endk,
            output_device_memory,'d');
        //                FindEdgeCutsCUDA(
        //                    input_device_memoryV,
        //                    device_scalesV,
        //                    device_scaleV,
        //                    w, h, d, m,
        //                    i, j, k, endi, endj, endk,
        //                    dual_edges_output,'d');
        CudaCheckReturn(cudaThreadSynchronize());
        //copy the block results back
        CudaCheckReturn(cudaMemcpy(dual_edges_output,
                                   output_device_memory,
                                   sizeof(Edge) *
                                   output_memory_size ,
                                   cudaMemcpyDeviceToHost));
        //call the axis edge kernel
        FindEdgeCutsCUDA<<<dimGrid3,dimBlock>>>(
            data,
            scales,
            scale,
            w, h, d, m,
            i, j, k, endi, endj, endk,
            output_device_memory,'a');
        //                FindEdgeCutsCUDA(
        //                    input_device_memoryV,
        //                    device_scalesV,
        //                    device_scaleV,
        //                    w, h, d, m,
        //                    i, j, k, endi, endj, endk,
        //                    axis_edges_output,'a');
        CudaCheckReturn(cudaThreadSynchronize());
        //copy the block results back
        CudaCheckReturn(cudaMemcpy(axis_edges_output,
                                   output_device_memory,
                                   sizeof(Edge) *
                                   output_memory_size ,
                                   cudaMemcpyDeviceToHost));
        CleaverCUDA::Edge* tmp;
        for (size_t ii = 0; ii < endi - i; ii++)
          for (size_t jj = 0; jj < endj - j; jj++)
            for (size_t kk = 0; kk < endk - k; kk++) {
              size_t cellIdx = i + ii + (j + jj)*wl + (k + kk)*wl*hl;
              size_t hostIdx = cellIdx + wl*hl;
              size_t deviceIdx = ii + jj * chunk + kk * chunk * chunk;
              //count the cuts
              for(size_t ee = 0; ee < 8; ee++) {
                //count diagonal (inner) edge cuts
                tmp = &inner_edges_output[deviceIdx*8+ee];
                if (cut_cells[cellIdx]) {
                  if ((tmp->isCut_eval & CleaverCUDA::kIsCut))
                    count++;
                } else tmp->isCut_eval = 0;
                if (ee < 3) {
                  //count dual edge cuts
                  tmp = &dual_edges_output[deviceIdx*3+ee];
                  bool include = cut_cells[cellIdx];
                  if (ee == 0 && cellIdx - 1 < max_cell)
                    include |= cut_cells[cellIdx - 1];
                  else if (ee == 1 && cellIdx - wl < max_cell)
                    include |= cut_cells[cellIdx - wl];
                  else if (ee  == 2 && cellIdx - wl*hl < max_cell)
                    include |= cut_cells[cellIdx - wl*hl];
                  //only include cut cells
                  if (!include)
                    tmp->isCut_eval = 0;
                  else if ((tmp->isCut_eval & CleaverCUDA::kIsCut))
                    count++;
                  //count axis edge cuts
                  tmp = &axis_edges_output[deviceIdx*3+ee];
                  include = cut_cells[cellIdx];
                  if (ee == 0) {
                    if(cellIdx - wl*hl < max_cell)
                      include |= cut_cells[cellIdx - wl*hl];
                    if(cellIdx - wl*hl - wl < max_cell)
                      include |= cut_cells[cellIdx - wl*hl -wl];
                    if(cellIdx - wl < max_cell)
                      include |= cut_cells[cellIdx - wl];
                  }
                  else if (ee == 1) {
                    if(cellIdx - wl*hl < max_cell)
                      include |= cut_cells[cellIdx - wl*hl];
                    if(cellIdx - wl*hl - 1 < max_cell)
                      include |= cut_cells[cellIdx - wl*hl - 1];
                    if(cellIdx - 1 < max_cell)
                      include |= cut_cells[cellIdx - 1];
                  }
                  else if (ee == 2) {
                    if(cellIdx - 1 < max_cell)
                      include |= cut_cells[cellIdx - 1];
                    if(cellIdx - 1 - wl < max_cell)
                      include |= cut_cells[cellIdx - 1 - wl];
                    if(cellIdx - wl < max_cell)
                      include |= cut_cells[cellIdx - wl];
                  }
                  //only include cut cells
                  if (!include)
                    tmp->isCut_eval = 0;
                  else if ((tmp->isCut_eval & CleaverCUDA::kIsCut))
                    count++;
                }
              }
              //inner edges copy to CPU
              memcpy(&inner_edges[hostIdx*8],
                     &inner_edges_output[deviceIdx*8], 8*sizeof(Edge));
              //dual edges copy to CPU
              memcpy(&dual_edges[hostIdx*3],
                     &dual_edges_output[deviceIdx*3], 3*sizeof(Edge));
              //axis edges copy to CPU
              memcpy(&axis_edges[hostIdx*3],
                     &axis_edges_output[deviceIdx*3], 3*sizeof(Edge));
            }
      }
    }
  }
  delete[] axis_edges_output;
  delete[] inner_edges_output;
  delete[] dual_edges_output;
  // Free GPU and reset
  CudaCheckReturn(cudaFree(data));
  CudaCheckReturn(cudaFree(scales));
  CudaCheckReturn(cudaFree(scale));
  CudaCheckReturn(cudaFree(output_device_memory));
  CudaCheckReturn(cudaDeviceReset());
  return count;
}
}
