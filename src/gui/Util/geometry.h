/**
 * @file geometry.h
 * This file contains the basic geometries and
 * geometry access for CleaverCUDA.
 * @version Feb 26, 2014
 * @author: Brig Bagley
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <array>
#include "cleaver_cuda.hh"

/** Macro to index into a 1D array with 3D dimensions. */
#define Idx3(a,b,c)           ((a) + w_ * (b) + w_ * h_ * (c))
/** Macro to index into a 1D array with 4D dimensions. */
#define Idx4(a,b,c,d)         \
    ((a) + w_ * (b) + w_ * h_ * (c) + w_ * h_ * d_ * (d))

/**
 * Simple Face struct used for CUDA faces.
 */
typedef struct SimpleFace {
  char isCut_eval = 0; //& 0xf0 for isCut, & 0x0f for eval
} SimpleFace;
/**
 * Simple tet struct used for CUDA tets.
 */
typedef struct SimpleTet {
  char isCut_eval = 0; //& 0xf0 for isCut, & 0x0f for eval
} SimpleTet;

/**
 * This is a list of numbers of geometries per cell.
 */
class Definitions {
 public:
  const static size_t kEdgesPerCell = 26;
  const static size_t kFacesPerCell = 36;
  const static size_t kTetsPerCell  = 24;

  enum tri_index {
    FUL, FUR, FUF, FUB,        //Inner upper  (0-3)
    FLL, FLR, FLF, FLB,        //Inner lower  (4-7)
    FFL, FFR, FBL, FBR,        //Inner column (8-11)
    FLUF, FLUB, FLLF, FLLB,    //Outer Left   (12-15)
    FRUF, FRUB, FRLF, FRLB,    //Outer Right  (16-19)
    FUFL, FUBL, FUFR, FUBR,    //Outer Upper  (20-23)
    FDFL, FDBL, FDFR, FDBR,    //Outer Down   (24-27)
    FFUL, FFUR, FFLL, FFLR,    //Outer Front  (28-31)
    FBUL, FBUR, FBLL, FBLR};   //Outer Back   (32-35)

  enum tet_index {
    TLU, TLL, TLF, TLB,        // Left  Face Tets (0-3)
    TRU, TRL, TRF, TRB,        // Right Face Tets (4-7)
    TFT, TFB, TFL, TFR,        // Front Face Tets (8-11)
    TBT, TBB, TBL, TBR,        // Back  Face Tets (12-15)
    TDF, TDB, TDL, TDR,        // Down  Face Tets (16-19)
    TUF, TUB, TUL, TUR};       // Upper Face Tets (20-23)

  enum cell_index {LEFT, RIGHT, UP, DOWN, FRONT, BACK, CENTER};
};
/**
 * SimpleGeometry holds pointers to all geometries and provides operations
 * to get the vertices, edges, and faces.
 */
class SimpleGeometry {
 public:
  /**
   * Constructor allocates all the memory needed for the geometry operated on.
   * @param w The width of the block.
   * @param h The hieght of the block.
   * @param d The depth of the block.
   * @param verbose Whether to print timing or not.
   */
  SimpleGeometry(size_t w, size_t h, size_t d,bool verbose);
  /**
   * Destructor deletes all of the memory declared.
   */
  ~SimpleGeometry();
  /**
   * Returns whether the pointers are all valid.
   * @return True if all pointers are valid.
   */
  bool Valid();
  /**
   * Gets the 3 point locations (vertices) of a face by it's number.
   * @param num The number of the face.
   * @param i The x location in the lattice.
   * @param j The y location in the lattice.
   * @param k The z location in the lattice.
   */
  std::array<std::array<float,3>,3> GetFaceVertices(
      Definitions::tri_index num,
      size_t i, size_t j, size_t k,
      std::array<float,3> scale);
  /**
   * Gets the 4 point locations (vertices) of a tet by it's number.
   * @param num The number of the tet.
   * @param i The x location in the lattice.
   * @param j The y location in the lattice.
   * @param k The z location in the lattice.
   */
  std::array<std::array<float,3>,4> GetTetVertices(
      Definitions::tet_index num,
      size_t i, size_t j, size_t k,
      std::array<float,3> scale);
  /**
   * Returns an array of 3 pointers to the 3 edges for the given face of a cell.
   * @param cell The number of the cell.
   * @param fnum  The number of the face.
   * @param tnum  The number of the tet.
   */
  std::array<CleaverCUDA::Edge*,3> GetFaceEdges(size_t cell,
                                         Definitions::tri_index fnum,
                                         Definitions::tet_index tnum);
  /**
   * Returns an array of 3 indices to the 3 edges for the given face of a cell.
   * @param num  The number of the face.
   */
  std::array<CleaverCUDA::edge_index,3> GetFaceEdgesNum(
      Definitions::tri_index num);
  /**
   * Returns an array of 4 pointers to the 4 faces for the given tet of a cell.
   * @param cell The number of the cell.
   * @param num  The number of the tet.
   */
  std::array<SimpleFace*,4> GetTetFaces(size_t cell,
                                        Definitions::tet_index num);
  /**
   * Returns an array of 6 pointers to the 6 edges for the given tet of a cell.
   * @param cell The number of the cell.
   * @param num  The number of the tet.
   */
  std::array<CleaverCUDA::Edge*,6> GetTetEdges(
      size_t cell, Definitions::tet_index n);
  /**
   * Returns an array of 4 indices of the 4 faces for the given tet of a cell.
   * @param num  The number of the tet.
   */
  std::array<Definitions::tri_index,4> GetTetFacesNum(
      Definitions::tet_index num);
  /**
   * Returns a pointer to the desired edge of a cell.
   * @param cell The cell number (<w*h*d).
   * @param edge The numbered edge (0-25) desired. 0-7 are the 8 inner
   * diagonals, 8-13 are the dual edges, and 14-25 are the face edges.
   * @return the pointer to the desired edge.
   */
  CleaverCUDA::Edge * GetEdge(int64_t cell,
                       CleaverCUDA::edge_index edge);
  /**
   * Returns the 3 pointers to the different edge arrays.
   * @return The 3 pointers to the different edge arrays.
   */
  std::array<CleaverCUDA::Edge *,3> GetEdgePointers();
  /**
   * Returns the sizes to the different edge arrays.
   * @return The sizes to the different edge arrays.
   */
  std::array<size_t,3> GetEdgePointersSize();
  /**
   * Returns a pointer to the desired face of a cell.
   * @param cell The cell number (<w*h*d).
   * @param face The face number.
   */
  SimpleFace * GetFace(int64_t cell,
                       Definitions::tri_index face);
  /**
   * Returns a pointer to the desired tet of a cell.
   * @param cell The cell number (<w*h*d).
   * @param tet The tet number.
   */
  SimpleTet* GetTet(int64_t cell,
                    Definitions::tet_index tet);
  /**
   * Tests all of the geometry to be sure it is correct.
   */
  void TestCells();
 private:
  /**
   * Creates a standard array from 3 size_t values.
   * @param i The x value.
   * @param j The y value.
   * @param k The z value.
   * @return The standard array of 3 size_t values.
   */
  static std::array<size_t,3> CreateArray(size_t i, size_t j, size_t k);
  /**
   * Returns the index for the inner edge of a cell from a cell index.
   * @param cell The cell index.
   * @param face The edge number to offset (0-7).
   * @return The index of the inner edge array we need for this cell and dir.
   */
  size_t GetInnerEdgeIdx(size_t cell, size_t edge);
  /**
   * Returns the location of the edge in the bisect edge array
   * given the cell of interest, & a direction (up, down, left, right,
   * forward, backward).
   * @param cell The cell number we wish to get a bisect edge from.
   * @param dir The direction of the set of 6 bisect edge to choose from.
   * @return The index of the bisect edge array we need for this cell and dir.
   */
  size_t GetDualEdgeIdx(size_t cell,
                        CleaverCUDA::edge_index edge);
  /**
   * Returns the location of the edge in the edge array
   * given the cell of interest, a direction (up, down, left, right,
   * forward, backward).
   * @param cell The cell number we wish to get a edge from.
   * @param dir The direction to get an adjacent edge. UP, RIGHT, and FRONT
   * are native to a cell. DOWN, LEFT, and BACK are adjacent cells.
   * @return The index of the cube edge array we need for this cell and dir.
   */
  size_t GetAxisEdgeIdx(size_t cell,
                        CleaverCUDA::edge_index edge);
  /**
   * Returns the index for the inner face of a cell from a cell index.
   * @param cell The cell index.
   * @param face The face number to offset (0-11).
   * @return The index of the inner face array we need for this cell and dir.
   */
  size_t GetInnerFaceIdx(size_t cell, size_t face);
  /**
   * Returns the location of the tet/ outer face in the tet/face array
   * given the cell of interest, a direction (up, down, left, right,
   * forward, backward), and a number (0-3).
   * @param cell The cell number we wish to get a tet/face from.
   * @param dir The direction of the set of 4 tet/face to choose from.
   * @param num The number of the tet/face (offset) from 0 to 4.
   * @return The index of the tet/ outer face array we need for
   * this cell and dir.
   */
  size_t GetOuterFaceIdx(size_t cell,
                         Definitions::cell_index dir,
                         unsigned char num);
  /**
   * Helper function that does the offsetting from a cell location to
   * the location of a vertex, edge, face, or tet in respective arrays.
   * This is only helpful for shared geometry in the 3 primal axes.
   * @param cell The cell number.
   * @param dir The direction we are looking at for the geometry in the cell.
   * @param multiple How many of the geometry is in each dimension of the cell.
   * @param num The number we are looking for in that direction (0-multiple).
   */
  size_t OffsetAlgorithm(size_t cell,
                         Definitions::cell_index dir,
                         unsigned char multiple,
                         unsigned char num);
  /** Pointer to the inner diagonal edges. */
  CleaverCUDA::Edge   * inner_edges_;
  /** Pointer to the dual edges. */
  CleaverCUDA::Edge   * dual_edges_ ;
  /** Pointer to the cube (axis) edges. */
  CleaverCUDA::Edge   * axis_edges_;
  /** Pointer to the inner faces. */
  SimpleFace   * inner_faces_;
  /** Pointer to the outer faces. */
  SimpleFace   * outer_faces_;
  /** Pointer to the tets. */
  SimpleTet    * tets_;
  /** The number of inner diagonal edges. */
  size_t    num_inner_edges_;
  /** The number of dual edges. */
  size_t    num_dual_edges_ ;
  /** The number of cube (axis) edges. */
  size_t    num_axis_edges_;
  /** The number of inner faces. */
  size_t    num_inner_faces_;
  /** The number of outer faces. */
  size_t    num_outer_faces_;
  /** The number of tets. */
  size_t    num_tets_;
  /** The lattice width */
  size_t w_;
  /** The lattice height */
  size_t h_;
  /** The lattice depth */
  size_t d_;
};

#endif /* GEOMETRY_H_ */
