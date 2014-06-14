/**
 * @file cleaver_cuda.hh
 * This file contains the GPU calls for Cleaver.
 * @version April 21, 2014
 * @author: Brig Bagley
 */
#ifndef CLEAVER_CUDA_H_
#define CLEAVER_CUDA_H_

namespace CleaverCUDA {
/**
 * Simple Edge struct used for CUDA Edges.
 */
typedef struct Edge {
  unsigned char isCut_eval;  //& 0xf0 for isCut, & 0x0f for eval
  float cut_loc[3];
} Edge;

enum edge_index {
  DULF, DULB, DURF, DURB,   //  ]_ Diagonal      (0-3)
  DLLF, DLLB, DLRF, DLRB,   //  ]  Edges         (4-7)
  CL, CR, CU, CD, CF, CB,   // Dual   Edges      (8-13)
  UL, UR, UF, UB,           // Top    Face Edges (14-17)
  LL, LR, LF, LB,           // Bottom Face Edges (18-21)
  FL, FR, BL, BR };         // Four Column Edges (22-25)

/** This constant is a mask for geometery evaluation */
static const unsigned char kIsEvaluated = 0x01;
/** This constant is a mask for an edge cut */
static const unsigned char kIsCut = 0x02;
/** This constant is a mask for a tet's stenciling */
static const unsigned char kIsStenciled = 0x04;
/** This constant is a mask for if a tet has > 2 cuts. */
static const unsigned char kHasStencil = 0x08;
/** This constant is a shift for an edge's material */
static const unsigned char kMaterial = 4;
/** The maximum number of materials allowed (as far as coloring) */
static const unsigned char kMaxMaterials = 6;


/** The constant size of a CUDA thread edge */
static const size_t kThreadSize = 8;
/** The constant size of a CUDA block edge */
static const size_t kBlockSize = 16;
/** The constant size of a lattice chunk edge */
static const size_t kChunkSize = kThreadSize*kBlockSize;

/**
 * The Kernel call for finding the dominant materials.
 * @param all_data The input data pointer.
 * @param scales The data scales for each material.
 * @param scale The user-defined resolution scale.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @param m The input data number of materials.
 * @param labels The pointer to the full array of data labels.
 * @param wl The label array width.
 * @param hl The label array height.
 * @param dl The label array depth.
 * @param The data, scales, and scale pointers on the device
 * as an array of 3 will be assigned at the end of the function.
 */
void CallCUDAMaxes(float *all_data,
                   float* scales, float* scale,
                   size_t w, size_t h, size_t d, size_t m,
                   char * labels,
                   size_t wl, size_t hl, size_t dl,
                   void* device_pointers[3]);

/**
 * The Kernel call for finding the edge cuts.
 * @param all_data The input data pointer.
 * @param scales The data scales for each material.
 * @param scale The user-defined resolution scale.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @param m The input data number of materials.
 * @param wl The label array width.
 * @param hl The label array height.
 * @param dl The label array depth.
 * @param inner_edges The pointer to the array of inner edges.
 * @param num_inner_edges The number of inner edges.
 * @param dual_edges The pointer to the array of dual edges.
 * @param num_dual_edges The number of dual edges.
 * @param axis_edges The pointer to the array of axis edges.
 * @param num_axis_edges The number of axis edges.
 * @return The number of edge cuts found.
 */
size_t CallCUDACuts(void* data,
                    void* scales, void* scale,
                    size_t w, size_t h, size_t d, size_t m,
                    size_t wl, size_t hl, size_t dl,
                    Edge* inner_edges,
                    Edge* dual_edges,
                    Edge* axis_edges,
                    bool* cut_cells);

/**
 * Finds the respective adjacent cell from a given one
 * that is useful when calculating values for edge endpoints.
 * @param num The edge of interest.
 * @param first Whether we want the first or second vertex on the edge.
 * @param i The x location in the block.
 * @param j The y location in the block.
 * @param k The z location in the block.
 * @param cell The array of the location of the relevant adjacent edge.
 */
void GetAdjacentCellFromEdgeCUDA(
    CleaverCUDA::edge_index num, bool first,
    size_t i, size_t j, size_t k, size_t* cell);

/**
 * Device function to get the cell center of a certain cell.
 * @param data_ The pointer to the material data.
 * @param scales The data specific scales.
 * @param scale The scaling of all of the data.
 * @param i The x location in the block.
 * @param j The y location in the block.
 * @param k The z location in the block.
 * @param m The material to grab.
 * @param num_mats The total number of materials.
 * @param find_max Whether to find the max material, or just the given one.
 * @param max_mat The char ptr to set the material number of the cell's center.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @return Returns the center vertex value.
 */
float GetCellCenterValueCUDA(
    float* data_, float* scales, float* scale,
    size_t i, size_t j, size_t k,
    size_t m, size_t num_mats,
    bool find_max, char* max_mat,
    size_t w, size_t h, size_t d);
/**
 * Indexes in to the material data to find a material value at a location.
 * This function is meant for GPU use.
 * @param data_ The pointer to the material data.
 * @param i The x location in the data.
 * @param j The y location in the data.
 * @param k The z location in the data.
 * @param mm The material to look at.
 * @param m the number of materials.
 * @param scales The data specific scales.
 * @param scale The scaling of all of the data.
 * @param ww_ The input data width.
 * @param hh_ The input data height.
 * @param dd_ The input data depth.
 *
 */
float DataTransformCUDA(float *data_,
                        float i, float j, float k,
                        size_t mm, size_t m, float* scales,
                        float* scale,
                        size_t ww_, size_t hh_, size_t dd_);

/**
 * Finds the material and value of that material at a given
 * vertex of an edge.
 * @param data_ The pointer to the material data.
 * @param scales The data specific scales.
 * @param scale The scaling of all of the data.
 * @param num The edge of interest.
 * @param first Whether we want the first or second vertex on the edge.
 * @param find_max Whether we want to find the max, or a particular material.
 * @param m The material to grab (must be in range) if applicable.
 * @param num_mats The total number of materials.
 * @param ret_mat The char pointer to set the material number of edge's end.
 * @param i The x location in the block.
 * @param j The y location in the block.
 * @param k The z location in the block.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @return Returns the edge's vertex material value.
 */
float GetEdgeMatAndValueAtEndpointCUDA(
    float* data_, float* scales, float* scale,
    CleaverCUDA::edge_index num, bool first, bool find_max,
    size_t m, size_t num_mats, char* ret_mat,
    size_t i, size_t j, size_t k,
    size_t w, size_t h, size_t d);

/**
 * Finds the material with the maximum value.
 * @param data_ The pointer to the material data.
 * @param scales The data specific scales.
 * @param scale The scaling of all of the data.
 * @param i The x location in the block.
 * @param j The y location in the block.
 * @param k The z location in the block.
 * @param num_mats The total number of materials.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @return The number of the material of highest value at this location.
 */
char GetLabelCUDA(float* data_, float* scales, float* scale,
              size_t i, size_t j, size_t k, char num_mats,
              size_t w, size_t h, size_t d);

/**
 * CUDA version of get edge vertices. Uses primitive arrays
 * (2 sets of 3 floats to set the vertex of interest.
 * @param num The number of the edge.
 * @param i The x location in the lattice.
 * @param j The x location in the lattice.
 * @param k The x location in the lattice.
 * @param scale The scaling of the point.
 * @param verts The pointer to the location to store the vertices.
 *
 */
void GetEdgeVerticesCUDA(
    CleaverCUDA::edge_index num,
    size_t i, size_t j, size_t k,
    float* scale, float verts[2][3]);

/**
 * Calculates the given edge.
 * @param data_ Pointer to the material data.
 * @param scales The pointer to the scales to use.
 * @param scale The pointer to the scale to use.
 * @param w The input data width.
 * @param h The input data height.
 * @param d The input data depth.
 * @param m The input data number of materials.
 * @param i The x location of the start of the chunk.
 * @param j The y location of the start of the chunk.
 * @param k The z location of the start of the chunk.
 * @param edge The pointer to the edge in memory.
 * @param num The number of the edge (0-25).
 */
void FindEdgeCutCUDA(
    float* data_,
    float* scales,
    float* scale,
    size_t w, size_t h, size_t d, size_t m,
    size_t i, size_t j, size_t k,
    Edge* edge, CleaverCUDA::edge_index num);

/**
 * Sets the 3 values of the cell array to the i j k values.
 * @param cell The pointer to the array of 3.
 * @param i The first value of the array to set.
 * @param j The second value of the array to set.
 * @param k The third value of the array to set.
 */
void SetArrayCUDA(size_t *cell, size_t i, size_t j, size_t k);

}
#endif /* CLEAVER_CUDA_H_ */
