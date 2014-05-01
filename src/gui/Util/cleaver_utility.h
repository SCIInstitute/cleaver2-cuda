/**
 * @file cleaver_utility.h
 * This file contains the main cleaver operations.
 * @version Mar. 28, 2014
 * @author: Brig Bagley
 */
#ifndef CLEAVER_UTILITY_H_
#define CLEAVER_UTILITY_H_
#include <cstring>
#include <vector>
#include "geometry.h"
#include "cleaver_cuda.hh"

class CleaverUtility {
 public:
  /**
   * Default constructor does nothing.
   */
  CleaverUtility();
  /**
   * Destructor cleans up all memory.
   */
  ~CleaverUtility();
  /**
   * Takes command line input and parses options for the program.
   * @param argc The number of command line arguments.
   * @param argv The list of arguments.
   */
  void ParseInput(int argc, char *argv[]);
  /**
   * Uses the NRRD library to read all of the files given into a
   * large array for manipulation and use later.
   * @return True on successful read, otherwise false.
   */
  bool LoadNRRDs();
  /**
   * Finds the cells that have a cut in them.
   */
  void FindCutCells();
  /**
   * Finds the materials with the greatest value for each cell.
   */
  void FindMaxes();
  /**
   * Gets the lattice width.
   * @return The lattice width.
   */
  size_t w();
  /**
   * Gets the lattice height.
   * @return The lattice width.
   */
  size_t h();
  /**
   * Gets the lattice depth.
   * @return The lattice width.
   */
  size_t d();
  /**
   * Gets the number of materials.
   * @return The number of materials.
   */
  size_t m();
  /**
   * Gets whether we want verbose info output.
   * @return True if we want verbose output.
   */
  bool verbose();
  /**
   * Runs through all of the edges of each cell to compute cuts.
   * @param geos Reference to the SimpleGeometery pointers.
   */
  void CalculateCuts(SimpleGeometry&  geos);
  /**
   * Runs through all of the faces of each cell to compute triples.
   * @param geos Reference to the SimpleGeometery pointers.
   */
  void CalculateTriples(SimpleGeometry&  geos);
  /**
   * Runs through all of the tets of each cell to compute quadruples.
   * @param geos Reference to the SimpleGeometery pointers.
   */
  void CalculateQuadruples(SimpleGeometry&  geos);
  /**
   * Takes all of the calculated data and creates stenciled faces.
   * @param geos Reference to the SimpleGeometery to set pointers to.
   */
  void StencilFaces(SimpleGeometry& geos);
  /**
   * Outputs the list of faces and verts to a PLY file for now.
   */
  void OutputToFile();
  /**
   * Completes all of the steps from loading NRRD to creating faces.
   * @param inputs The file inputs to load the NRRDs.
   * @param scales The scaled resolution factors.
   * @param abs The absolute resolution dimensions.
   * @param useScale Whether to use the scale res. or not.
   * @param useAbs Whether to use the absolute res. or not.
   * @param useGPU Whether to use the GPU or not.
   * @return The pair of the vector of vertices and vector of faces.
   */
  std::pair<std::vector<std::array<float,3>>,
  std::vector<std::array<size_t,4>>> GetVertsFacesFromNRRD(
      std::vector<std::string> &files,
      float scales[3],
      std::array<size_t,3> &abs,
      bool useScale, bool useAbs, bool useGPU, bool addAir);
  /**
   * Parses command line input, but does the same as overloaded
   * method above.
   * @param argc The number of arguments.
   * @param argv The list of arguments.
   */
  void GetVertsFacesFromNRRD(int argc, char* argv[]);
 private:
  /**
   * Adds a new face to the mesh. If a vertex already exists in the
   * list, it isn't added again. Otherwise, it is, and the correct
   * vertex numbers are given for the faces.
   * @param face The face (3 verts) to add.
   * @param mat The material of the face
   */
  void AddFace(std::array<std::array<float,3>,3> face, char mat);
  /**
   * Finds all of the transitional faces in a tet.
   * @param tet The pointer to the tet.
   * @param num The tet number.
   * @param i The x location in the block.
   * @param j The y location in the block.
   * @param k The z location in the block.
   * @param geos The geometry pointers.
   */
  void AccumulateVertsAndFaces(
      SimpleTet * tet, Definitions::tet_index num,
      size_t i, size_t j, size_t k, SimpleGeometry& geos);
  /**
   * Computes the triple of a face.
   * @param face The edge we are working on.
   * @param num The face index.
   * @param i The x location in the block.
   * @param j The y location in the block.
   * @param k The z location in the block.
   * @param geos The pointers to all the geometry.
   */
  void CalculateTriple(SimpleFace *face, Definitions::tri_index num,
                       size_t i, size_t j, size_t k, SimpleGeometry& geos);
  /**
   * Computes the quadruple of a tet.
   * @param tet The edge we are working on.
   * @param num The tet index.
   * @param i The x location in the block.
   * @param j The y location in the block.
   * @param k The z location in the block.
   * @param geos The pointers to all the geometry.
   */
  void CalculateQuadruple(SimpleTet *tet, Definitions::tet_index num,
                       size_t i, size_t j, size_t k, SimpleGeometry& geos);
  /**
   * Determines where a triple point is on a face.
   * @param points The three points in space.
   * @return Returns a 3 float array of the triple point in x,y,z.
   */
  std::array<float,3> FindTriplePoint(
      std::array<std::array<float,3>,3> points);
  /**
   * Prints information on how to use the program.
   */
  void PrintUsage();
  /**
   * Prints information about the program and how to use it.
   */
  void PrintHelp();
  /** The list of input files. */
  std::vector<std::string> inputs_;
  /** The name of the ouput file. */
  std::string output_;
  /** The name of the format to output as. */
  std::string format_;
  /** Whether to have program info printed during execution. */
  bool verbose_;
  /** Whether to use absolute resolution or not */
  bool absolute_resolution_;
  /** Whether to use scaled resolution or not */
  bool scaled_resolution_;
  /** The 3 absolute resolution values. */
  std::array<size_t,3> res_;
  /** The 3 scaled resolution values. */
  float scale_[3];
  /** The pointer to all of the float data. */
  float * data_;
  /** The pointer to the data scales.*/
  std::vector<std::array<float,3>> scales_;
  /** The lattice width. */
  size_t w_;
  /** The lattice height. */
  size_t h_;
  /** The lattice depth. */
  size_t d_;
  /** The data width. */
  size_t ww_;
  /** The data height. */
  size_t hh_;
  /** The data depth. */
  size_t dd_;
  /** The number of materials. */
  size_t m_;
  /** An array of labels for which material is max at each cell. */
  char * labels_;
  /** An array of bools for which cells contain cuts. */
  bool * cut_cells_;
  /** The list of vertices accumulated after stenciling. */
  std::vector<std::array<float,3>> verts_;
  /** The list of faces accumulated after stenciling. */
  std::vector<std::array<size_t,4>> faces_;
  /** Whether to call the GPU or not. */
  bool use_GPU_;
  /** An array of 3 pointers used more than once for GPU */
  void* device_pointers_[3];
  /** A pointer to the data scales for GPU use */
  float* scalesP_;
  /** A bool to decide to add "air" material or not. */
  bool addAir_;

};

#endif /* CLEAVER_UTILITY_H_ */
