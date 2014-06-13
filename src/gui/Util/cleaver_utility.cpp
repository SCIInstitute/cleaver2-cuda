/**
 * @file cleaver_utility.cpp
 * This file contains the main cleaver operations.
 * @version Mar. 28, 2014
 * @author: Brig Bagley
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstring>
#include <teem/nrrd.h>
#include "cleaver_utility.h"
#ifdef CUDA_FOUND
#include "cleaver_cuda.hh"
#endif

CleaverUtility::CleaverUtility()
: output_("output"),
  format_("tetgen"),
  verbose_(true),
  absolute_resolution_(false),
  scaled_resolution_(false),
  data_(nullptr),
  w_(0),
  h_(0),
  d_(0),
  ww_(0),
  hh_(0),
  dd_(0),
  m_(0),
  labels_(nullptr),
  cut_cells_(nullptr),
  use_GPU_(true),
  addAir_(false) {
  device_pointers_[0] = device_pointers_[1] = device_pointers_[2] = NULL;
  scale_[0] = scale_[1] = scale_[2] = 1.;
  scalesP_ = nullptr;
}

CleaverUtility::~CleaverUtility() {
  if (verbose_)
    std::cout << "Clean up..." << std::endl;
  delete[] data_;
  delete[] labels_;
  delete[] cut_cells_;
  delete[] scalesP_;
  if (verbose_)
    std::cout << "Done." << std::endl;
}

void CleaverUtility::ParseInput(int argc, char *argv[])
{
  if(argc < 2) PrintUsage();
  int a = 0;
  std::string token;
  while(++a < argc) {
    token = argv[a];
    //------------------------
    //  Parse Help Flag
    //------------------------
    if(token.compare("-h") == 0)
      PrintHelp();
    //------------------------
    //  Parse GPU Flag
    //------------------------
    else if(token.compare("--force-CPU") == 0)
      use_GPU_ = false;
    //------------------------
    //  Parse Air Flag
    //------------------------
    else if(token.compare("--add-air") == 0)
      addAir_ = true;
    //------------------------
    //  Parse Silent Flag
    //------------------------
    else if(token.compare("-s") == 0)
      verbose_ = false;
    //--------------------------
    //   Parse Input Flag
    //--------------------------
    else if(token.compare("-i") == 0) {
      while(++a < argc) {
        token = argv[a];
        if (token.find(".nrrd") == std::string::npos){
          token = argv[--a];
          break;
        }
        inputs_.push_back(token);
      }
      m_ = inputs_.size();
    }
    //--------------------------
    //   Parse Output Flag
    //--------------------------
    else if(token.compare("-o") == 0) {
      if (++a >= argc) {
        std::cerr << "No output provided with '-o'." << std::endl;
        PrintUsage();
      }
      output_ = argv[a];
    }
    //--------------------------
    //   Parse Format Flag
    //--------------------------
    else if(token.compare("-f") == 0) {
      if (++a >= argc) {
        std::cerr << "No format provided with '-f'." << std::endl;
        PrintUsage();
      }
      format_ = argv[a];
      if(!(format_.compare("scirun") == 0 ||
          format_.compare("tetgen") == 0  ||
          format_.compare("matlab") == 0  )) {
        std::cerr << "invalid output format (-f): " << format_ << std::endl;
        PrintUsage();
      }
    }
    //-------------------------
    //   Absolute Resolution Flag
    //-------------------------
    else if(token.compare("-ra") == 0) {
      absolute_resolution_ = true;
      if(a+3 < argc) {
        res_[0] = atoi(argv[++a]);
        res_[1] = atoi(argv[++a]);
        res_[2] = atoi(argv[++a]);
        if (res_[0] == 0 || res_[1] == 0 || res_[2] == 0) {
          std::cerr << "Invalid -ra parameters." << std::endl;
          PrintUsage();
        }
      } else {
        std::cerr << "Invalid -ra parameters." << std::endl;
        PrintUsage();
      }
    }
    //-------------------------
    //   Scaled Resolution Flag
    //-------------------------
    else if(token.compare("-rs") == 0) {
      scaled_resolution_ = true;
      if(a+3 < argc) {
        scale_[0] = atof(argv[++a]);
        scale_[1] = atof(argv[++a]);
        scale_[2] = atof(argv[++a]);
        if (scale_[0] == 0 || scale_[1] == 0 || scale_[2] == 0) {
          std::cerr << "Invalid -ra parameters." << std::endl;
          PrintUsage();
        }
      } else {
        std::cerr << "Invalid -rs parameters." << std::endl;
        PrintUsage();
      }
    } else {
      std::cerr << "Unknown option flag: " << token << std::endl;
      PrintUsage();
    }
  }
  if (inputs_.size() == 0) {
    std::cerr << "No inputs were given." << std::endl;
    PrintUsage();
  }
}

void CleaverUtility::PrintUsage()
{
  std::cerr << "usage: cleaver -i INPUT1.nrrd [INPUT2.nrrd ...] " <<
      "[-h] [-s] [-p] [-as VAL] [-al VAL]  [-o NAME] [-f FORMAT] " <<
      "[-ra X Y Z] [-rs X Y Z] \n\n" << std::endl;
  std::cerr << "   -i   input filenames   minimimum=1   " << std::endl;
  std::cerr << "   -h   help                            " << std::endl;
  std::cerr << "   -s   silent mode       default=off   " << std::endl;
  //  std::cerr << "   -as  alpha short       default=0.357 " << std::endl;
  //  std::cerr << "   -al  alpha long        default=0.203 " << std::endl;
  std::cerr << "   -o   output filename   default=output" << std::endl;
  //  std::cerr << "   -f   output format     default=tetgen" << std::endl;
  std::cerr << "   -ra  x y z             absolute resolution" << std::endl;
  std::cerr << "   -rs  x y z             scaled resolution" << std::endl;
  std::cerr << "   --force-CPU            forces CPU only" << std::endl;
  std::cerr << "   --add-air              Adds \"Air\" Material" << std::endl;

  std::cerr << std::endl;

  //  std::cerr << "   Valid Parameters:                " << std::endl;
  //  std::cerr << "        alpha short       0.0 to 1.0" << std::endl;
  //  std::cerr << "        alpha long        0.0 to 1.0" << std::endl;
  //  std::cerr << "        tetmesh formats   tetgen, scirun, matlab"
  //      << std::endl << std::endl;

  std::cerr << "Examples:" << std::endl;
  std::cerr << "cleaver -h                                         "
      << "print help guide" << std::endl;
  std::cerr << "cleaver -i mat1.nrrd mat2.nrrd                     "
      << "basic use" << std::endl;
  std::cerr << "cleaver -i mat1.nrrd mat2.nrrd -rs .5 .5 .5        "
      << "scale resolution" << std::endl;
  std::cerr << "cleaver -i mat1.nrrd mat2.nrrd -o mesh             "
      << "specify target name" << std::endl;
  std::cerr << "cleaver -i mat1.nrrd mat2.nrrd -ra 100 100 80      "
      << "absolute resolution" << std::endl;
  std::cerr << "cleaver -s -i mat1.nrrd mat2.nrrd                  "
      << "silent mode" << std::endl << std::endl;

  exit(0);
}

void CleaverUtility::PrintHelp()
{
  /********Update Version info Here ****************/
  std::string VersionNumber = "1.5.4";
  std::string VersionDate = "Dec 20, 2013";
  std::string Version = std::string("Version") + " " +
      VersionNumber + " " + VersionDate;

  std::cerr << "Cleaver" << std::endl;
  std::cerr << "A Conforming, Multimaterial, Tetrahedral Mesh Generator" << std::endl;
  std::cerr << "with Bounded Element Quality." << std::endl;
  std::cerr << Version << std::endl << std::endl;

  std::cerr << "Copyright (C) 2013" << std::endl;
  std::cerr << "Jonathan Bronson" << std::endl;
  std::cerr << "bronson@sci.utah.edu" << std::endl;
  std::cerr << "Scientific Computing & Imaging Institute" << std::endl;
  std::cerr << "University of Utah, School of Computing" << std::endl;
  std::cerr << "Salt Lake City, Utah" << std::endl << std::endl;

  std::cerr << "What Can Cleaver Do?" << std::endl << std::endl;
  std::cerr << "  " << "Cleaver generates conforming tetrahedral meshes for" << std::endl;
  std::cerr << "  " << "multimaterial or multiphase volumetric data. Both  " << std::endl;
  std::cerr << "  " << "geometric accuracy and element quality are bounded." << std::endl;
  std::cerr << "  " << "The method is a stencil-based approach, and relies " << std::endl;
  std::cerr << "  " << "on an octree structure to provide a coarse level of" << std::endl;
  std::cerr << "  " << "grading in regions of homogeneity.                 " << std::endl;

  std::cerr << "What does Cleaver use as input?" << std::endl << std::endl;
  std::cerr << "  " << "The cleaving algorithm works by utilizing indicator" << std::endl;
  std::cerr << "  " << "functions. These functions indicate the strength or" << std::endl;
  std::cerr << "  " << "relative presence of a particular material. At each" << std::endl;
  std::cerr << "  " << "point, only the material with the largest indicator" << std::endl;
  std::cerr << "  " << "value is considered present. In practice, inside-  " << std::endl;
  std::cerr << "  " << "outside and distance functions are most common.    " << std::endl << std::endl;

  std::cerr << "What is the input format?" << std::endl;
  std::cerr << "  " << "Cleaver takes the Nearly Raw Raster Data format, or" << std::endl;
  std::cerr << "  " << "NRRD, as input. Information on the format is avail-" << std::endl;
  std::cerr << "  " << "able at the Teem website:                          " << std::endl;
  std::cerr << "  " << "http://teem.sourceforge.net/nrrd/format.html       " << std::endl << std::endl;

  std::cerr << "What is the output format?" << std::endl;
  std::cerr << "  " << "Cleaver can output into three mesh formats." << std::endl;
  std::cerr << "  " << "1) TetGen:  .node, .ele                    " << std::endl;
  std::cerr << "  " << "2) SCIRun:  .pts,  .elem, .txt             " << std::endl;
  std::cerr << "  " << "3) Matlab:  .mat                           " << std::endl;
  std::cerr << "  " <<  std::endl;
  std::cerr << "  " << "In addition, Cleaver outputs a .info file  " << std::endl;
  std::cerr << "  " << "with more details about the output mesh.   " << std::endl << std::endl;

  PrintUsage();

  std::cerr << std::endl << std::endl;
}

bool CleaverUtility::LoadNRRDs() {
  clock_t start = clock();
  m_ = inputs_.size();
  if (m_ <= 0) return false;
  scales_.clear();
  //--------------------------------------------
  //  Read only headers of each file
  //-------------------------------------------
  std::vector<Nrrd*> nins;
  for(unsigned int i=0; i < inputs_.size(); i++) {
    // create empty nrrd file container
    nins.push_back(nrrdNew());

    // load only the header, not the data
    NrrdIoState *nio = nrrdIoStateNew();
    nrrdIoStateSet(nio, nrrdIoStateSkipData, AIR_TRUE);

    // Read in the Header
    if(verbose_) std::cout <<
        "Reading File: " << inputs_[i] << std::endl;
    if(nrrdLoad(nins[i], inputs_[i].c_str(), nio))
    {
      char *err = biffGetDone(NRRD);
      std::cerr << "Trouble Reading File: " <<
          inputs_[i] << " : " << err << std::endl;
      free(err);
      nio = nrrdIoStateNix(nio);
      continue;
    }

    // Done with nrrdIoState
    nio = nrrdIoStateNix(nio);
    delete nio;

    if(nins[i]->dim != 3)
    {
      std::cerr << "Fatal Error: volume dimension " <<
          nins[i]->dim << ", expected 3." << std::endl;
      for(size_t j = 0; j < i; j++)
        nrrdNuke(nins[j]);
      return false;
    }
  }
  //-----------------------------------
  // Verify contents match
  //-----------------------------------
  bool match = true;

  for(unsigned i=1; i < inputs_.size(); i++)
  {
    if(nins[i]->dim != nins[0]->dim){
      std::cerr << "Error, file " << i <<
          " # dims don't match." << std::endl;
      match = false;
    }

    for(unsigned int j=0; j < nins[0]->dim-1; j++)
      if(nins[i]->axis[j].size != nins[0]->axis[j].size)
      {
        std::cerr << "Error, file " << j <<
            " dimensions don't match." << std::endl;
        match = false;
        break;
      }

    if((int)nrrdElementNumber(nins[i]) != (int)nrrdElementNumber(nins[0]))
    {
      std::cerr << "Error, file " << i <<
          " does not contain expected number of elements" << std::endl;
      std::cerr << "Expected: " <<
          (int)nrrdElementNumber(nins[0]) << std::endl;
      std::cerr << "Found:    " <<
          (int)nrrdElementNumber(nins[i]) << std::endl;

      match = false;
      break;
    }
  }

  if(!match)
  {
    for (size_t i = 0; i < nins.size(); i++)
      nrrdNuke(nins[i]);
    return false;
  }


  //-----------------------------------
  // Save the dimensions
  //-----------------------------------
  w_ = ww_ = nins[0]->axis[0].size;
  h_ = hh_ = nins[0]->axis[1].size;
  d_ = dd_ = nins[0]->axis[2].size;

  if(verbose_)
    std::cout << "Input Dimensions: " << w_ <<
    " x " << h_ << " x " << d_ << std::endl;

  //---------------------------------------
  //     Allocate Sufficient Memory
  //---------------------------------------
  size_t total_cells = ww_*hh_*dd_;
  delete[] data_;
  data_ = new float[total_cells*m_];
  if (!data_) {
    std::cerr << "failed to allocate data memory." << std::endl;
    return false;
  }
  //--------------------------------------
  //  Deferred  Data Load/Copy
  //----------------------------------------
  for(unsigned int f=0; f < inputs_.size(); f++)
  {
    // load nrrd data, reading memory into data array
    if(nrrdLoad(nins[f], inputs_[f].c_str(), NULL))
    {
      char *err = biffGetDone(NRRD);
      std::cerr << "trouble reading data in file: " <<
          inputs_[f] << " : " << err << std::endl;
      free(err);
      for (size_t i = 0; i < nins.size(); i++)
        nrrdNuke(nins[i]);
      return false;
    }
    float (*lup)(const void *, size_t I);
    lup = nrrdFLookup[nins[f]->type];

    // cast and copy into large array
    int s=0;
    for(size_t k=0; k < d_; k++){
      for(size_t j=0; j < h_; j++){
        for(size_t i=0; i < w_; i++){
          data_[i + j*w_ + k*w_*h_ +
                f*total_cells ] = lup(nins[f]->data, s++);
        }
      }
    }

    // set scale
    float xs = ((Nrrd*)nins[f])->axis[0].spacing;
    float ys = ((Nrrd*)nins[f])->axis[1].spacing;
    float zs = ((Nrrd*)nins[f])->axis[2].spacing;

    // handle NaN cases
    if(xs != xs) xs = 1.;
    if(ys != ys) ys = 1.;
    if(zs != zs) zs = 1.;
    std::array<float,3> scale = {{ xs, ys, zs }};
    scales_.push_back(scale);
  }

  if (addAir_ || inputs_.size() == 1) {
    std::cerr << "Attempting transition mesh with 1 material..."
        << std::endl;
    inputs_.resize(inputs_.size()+1);
    m_ += 1;
    scales_.resize(scales_.size()+1);
    scales_[scales_.size()-1] = scales_[0];
    float* data2 = new float[total_cells*m_];
    for(size_t k=0; k < d_; k++){
      for(size_t j=0; j < h_; j++){
        for(size_t i=0; i < w_; i++){
          for(size_t h=0; h < m_; h++){
            size_t base = i + j*w_ + k*w_*h_;
            size_t base2 = base + h*total_cells;
            if (h < m_ - 1)
              data2[base2] = data_[base2];
            else {
              float mx = data_[base];
              for(size_t t = 1; t < m_ - 1; t++)
                mx = std::max(mx,data_[base + t*total_cells]);
              data2[base2] = -mx;
            }
          }
        }
      }
    }
    delete[] data_;
    data_ = data2;
  }




  delete[] scalesP_;
  scalesP_ = new float[scales_.size()*3];
  for (size_t u = 0; u < scales_.size(); u++)
    for (size_t v = 0; v < 3; v++)
      scalesP_[u*3 + v] = scales_.at(u).at(v);
  //set absolute resolution
  if(absolute_resolution_ && !scaled_resolution_) {
    w_ = res_[0];
    h_ = res_[1];
    d_ = res_[2];
    scale_[0] = static_cast<float>(ww_) / static_cast<float>(w_);
    scale_[1] = static_cast<float>(hh_) / static_cast<float>(h_);
    scale_[2] = static_cast<float>(dd_) / static_cast<float>(d_);
    scaled_resolution_ = true;
  }
  //scale the resolution
  else if(scaled_resolution_ && ! absolute_resolution_) {
    w_ *= scale_[0];
    h_ *= scale_[1];
    d_ *= scale_[2];
    scale_[0] = 1. / scale_[0];
    scale_[1] = 1. / scale_[1];
    scale_[2] = 1. / scale_[2];
  }
  //scale the resolution and absolute resolution
  else if(scaled_resolution_ && absolute_resolution_) {
    w_ = scale_[0] * res_[0];
    h_ = scale_[1] * res_[1];
    d_ = scale_[2] * res_[2];
    scale_[0] = static_cast<float>(ww_) / static_cast<float>(w_);
    scale_[1] = static_cast<float>(hh_) / static_cast<float>(h_);
    scale_[2] = static_cast<float>(dd_) / static_cast<float>(d_);
  }
  for (size_t i = 0; i < nins.size(); i++)
    nrrdNuke(nins[i]);
  //always pad.
  m_ = inputs_.size()+1;
  w_++;
  h_++;
  d_++;
  //  Allocate memory to find max material in each cell
  delete[] labels_;
  labels_ = new char[w_*h_*d_];
  if (!labels_) {
    std::cerr << "Failed to allocate labels char array. " << std::endl;
    return false;
  }
  //  Set up a list of booleans for cut/buffer cells.
  delete[] cut_cells_;
  cut_cells_ = new bool[w_*h_*d_];
  if (!cut_cells_) {
    std::cerr << "Failed to allocate cuts bool array. " << std::endl;
    return false;
  }
  //set defaults to false
  for (size_t i = 0; i < w_*h_*d_; i++)
    cut_cells_[i] = false;
  double duration = ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "\nRead Data:\t\t\t" << duration << " sec." << std::endl;
  return true;
}

void CleaverUtility::FindMaxes() {
  clock_t start = clock();
#ifdef CUDA_FOUND
  if (use_GPU_) {
    CleaverCUDA::CallCUDAMaxes(data_,scalesP_,scale_,
                               ww_,hh_,dd_,m_,labels_,
                               w_,h_,d_,device_pointers_);
    double duration = ((double)clock() -
        (double)start) / (double)CLOCKS_PER_SEC;
    if (verbose_)
      std::cout << "Found maxes:\t\t\t" <<
      duration << " sec." << std::endl;
    return;
  }
#endif
  for (size_t i = 0; i < w_; i++)
    for (size_t j = 0; j < h_; j++)
      for (size_t k = 0; k < d_; k++) {
        int max = 0;
        float max_val = CleaverCUDA::DataTransformCUDA(
            data_,static_cast<float>(i),
            static_cast<float>(j),
            static_cast<float>(k),
            0,m_,scalesP_,scale_,ww_,hh_,dd_);
        for (size_t l = 1; l < m_; l++) {
          float tmp;
          if ((tmp = CleaverCUDA::DataTransformCUDA(
              data_,static_cast<float>(i),
              static_cast<float>(j),
              static_cast<float>(k),
              l,m_,scalesP_,scale_,ww_,hh_,dd_)) > max_val) {
            max_val = tmp;
            max = l;
          }
        }
        labels_[Idx3(i,j,k)] = max;
      }
  double duration = ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "Found maxes:\t\t\t" << duration << " sec." << std::endl;
}

void CleaverUtility::FindCutCells() {
  clock_t start = clock();
  for (size_t i = 0; i < w_; i++)
    for (size_t j = 0; j < h_; j++)
      for (size_t k = 0; k < d_; k++) {
        if ((i+1 >= w_) || (j+1 >= h_) || (k+1 >= d_)) continue;
        if( labels_[Idx3(i,j,k)]     != labels_[Idx3(i+1,j,k)]     ||
            labels_[Idx3(i,j,k+1)]   != labels_[Idx3(i+1,j,k+1)]   ||
            labels_[Idx3(i,j,k)]     != labels_[Idx3(i,j,k+1)]     ||
            labels_[Idx3(i+1,j,k)]   != labels_[Idx3(i+1,j,k+1)]   ||

            labels_[Idx3(i,j+1,k)]   != labels_[Idx3(i+1,j+1,k)]   ||
            labels_[Idx3(i,j+1,k+1)] != labels_[Idx3(i+1,j+1,k+1)] ||
            labels_[Idx3(i,j+1,k)]   != labels_[Idx3(i,j+1,k+1)]   ||
            labels_[Idx3(i+1,j+1,k)] != labels_[Idx3(i+1,j+1,k+1)] ||

            labels_[Idx3(i,j,k)]     != labels_[Idx3(i,j+1,k)]     ||
            labels_[Idx3(i+1,j,k)]   != labels_[Idx3(i+1,j+1,k)]   ||
            labels_[Idx3(i,j,k+1)]   != labels_[Idx3(i,j+1,k+1)]   ||
            labels_[Idx3(i+1,j,k+1)] != labels_[Idx3(i+1,j+1,k+1)]) {
          cut_cells_[Idx3(i,j,k)] = true;
        }
      }
  //count cells
  size_t count_cuts = 0;
  for (size_t i = 0; i < w_*h_*d_; i++)
    if (cut_cells_[i]) count_cuts++;
  double duration = ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "Found cut cells: (" << count_cuts << ")\t" <<
    ((count_cuts < 10000)?"\t":"") <<
    duration << " sec." << std::endl;
}

size_t CleaverUtility::w() { return w_; }

size_t CleaverUtility::h() { return h_; }

size_t CleaverUtility::d() { return d_; }

size_t CleaverUtility::m() { return m_; }

bool CleaverUtility::verbose() { return verbose_; }

void CleaverUtility::CalculateCuts(SimpleGeometry&  geos) {
  clock_t start = clock();
#ifdef CUDA_FOUND
  if (use_GPU_) {
    size_t count_cuts =
        CleaverCUDA::CallCUDACuts(device_pointers_[0],
                                  device_pointers_[1],
                                  device_pointers_[2],
                                  ww_,hh_,dd_,m_,
                                  w_,h_,d_,
                                  geos.GetEdgePointers()[0],
                                  geos.GetEdgePointers()[1],
                                  geos.GetEdgePointers()[2],
                                  cut_cells_);
    //            CleaverCUDA::CallCUDACuts(data_,
    //                                      scalesP_,
    //                                      scale_,
    //                                      ww_,hh_,dd_,m_,
    //                                      w_,h_,d_,
    //                                      geos.GetEdgePointers()[0],
    //                                      geos.GetEdgePointers()[1],
    //                                      geos.GetEdgePointers()[2],
    //                                      cut_cells_);
    //    for (size_t tt = 0; tt < 3; tt++)
    //      for (size_t t = 0; t < geos.GetEdgePointersSize()[tt]; t++)
    //        if((geos.GetEdgePointers()[tt][t].isCut_eval & CleaverCUDA::kIsCut))
    //          std::cout << tt << " - " <<t << " : " <<
    //          (int)geos.GetEdgePointers()[tt][t].isCut_eval
    //          << " -- " << geos.GetEdgePointers()[tt][t].cut_loc[0] << ", "
    //          << geos.GetEdgePointers()[tt][t].cut_loc[1] << ", " <<
    //          geos.GetEdgePointers()[tt][t].cut_loc[2] << "\n";
    double duration = ((double)clock() - (double)start) /
        (double)CLOCKS_PER_SEC;
    if (verbose_)
      std::cout << "Found all Cuts: (" << count_cuts << ")\t" <<
      ((count_cuts < 100000)?"\t":"") << duration << " sec." << std::endl;
    return;
  }
#endif
  size_t count_cuts = 0;
  for(size_t i = 0; i < w_; i ++)
    for(size_t j = 0; j < h_; j ++)
      for(size_t k = 0; k < d_; k ++) {
        if(cut_cells_[Idx3(i,j,k)]) {
          for(size_t l = 0; l < Definitions::kEdgesPerCell; l++) {
            CleaverCUDA::Edge* edge =
                geos.GetEdge(Idx3(i,j,k),(CleaverCUDA::edge_index)l);
            if(!(edge->isCut_eval & CleaverCUDA::kIsEvaluated)) {
              CleaverCUDA::FindEdgeCutCUDA(
                  data_,scalesP_,scale_,ww_,hh_,dd_,m_,i,j,k,edge,
                  (CleaverCUDA::edge_index)l);
              if ((edge->isCut_eval & CleaverCUDA::kIsCut)) count_cuts++;
            }
          }
        }
      }
  //  for (size_t tt = 0; tt < 3; tt++)
  //    for (size_t t = 0; t < geos.GetEdgePointersSize()[tt]; t++)
  //      if((geos.GetEdgePointers()[tt][t].isCut_eval & CleaverCUDA::kIsCut))
  //        std::cout << tt << " - " <<t << " : " <<
  //        (int)geos.GetEdgePointers()[tt][t].isCut_eval
  //        << " -- " << geos.GetEdgePointers()[tt][t].cut_loc[0] << ", "
  //        << geos.GetEdgePointers()[tt][t].cut_loc[1] << ", " <<
  //        geos.GetEdgePointers()[tt][t].cut_loc[2] << "\n";
  double duration = ((double)clock() -
      (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "Found all Cuts: (" << count_cuts << ")\t" <<
    ((count_cuts < 100000)?"\t":"") << duration << " sec." << std::endl;
}

std::array<float,3> CleaverUtility::FindTriplePoint(
    std::array<std::array<float,3>,3> points) {
  //  float epsilon = 1e-5;
  //  float offset = 1e-2;
  //  //get the values of all materials at each point.
  //  std::array<std::vector<float>,3> mats;
  //  for (auto& a : mats) a.resize(m);
  //  for (size_t i = 0; i < 3; i++)
  //    for(size_t j = 0; j < m; j++)
  //      mats[i][j] = data_transform(data,points[i][0], points[i][1],
  //                                  points[i][2], m,w,h,d);
  //  //get the value of all materials at a point very close to point 1
  //  // on the line between point 1 and point 2. (X)
  //  float off_pt[3];
  //  for (auto a = 0; a < 3; a++)
  //    off_pt[a] = offset * points[0][a] + (1.-offset) * points[1][a];
  //  std::vector<float> off_vals;
  //  for (auto a = 0; a < m; a++)
  //    off_vals.push_back(data_transform(data,off_pt[0], off_pt[1],
  //                                      off_pt[2], m,w,h,d));
  //  //find the points on AB and AX where materials are equal. TODO
  //
  //  //repeat the steps above with AC and CB
  //
  //  //intersect the 2 lines from both sets of points (the answer.)
  //
  //  //check that the material values are all very close at this point.
  //  return std::array<float,3>();
  std::array<float,3> ans;
  for (auto a = 0; a < 3; a++)
    ans[a] = (points[0][a] + points[1][a] + points[2][a]) / 3.;
  return ans;
}

void CleaverUtility::CalculateTriples(SimpleGeometry& geos) {
  clock_t start = clock();
  size_t count_trips = 0;
  for(size_t i = 0; i < w_; i ++)
    for(size_t j = 0; j < h_; j ++)
      for(size_t k = 0; k < d_; k ++) {
        if(cut_cells_[Idx3(i,j,k)]) {
          for(size_t l = 0; l < Definitions::kFacesPerCell; l++) {
            SimpleFace* face = geos.GetFace(Idx3(i,j,k),
                                            (Definitions::tri_index)l);
            if(!(face->isCut_eval & CleaverCUDA::kIsEvaluated)) {
              CalculateTriple(face,(Definitions::tri_index)l,i,j,k,geos);
              if ((face->isCut_eval & CleaverCUDA::kIsCut)) count_trips++;
            }
          }
        }
      }
  double duration = ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "Found all Triples: (" << count_trips << ")\t" <<
    duration << " sec." << std::endl;
}

void CleaverUtility::CalculateTriple(SimpleFace *face,
                                     Definitions::tri_index num,
                                     size_t i, size_t j, size_t k,
                                     SimpleGeometry& geos)  {
  face->isCut_eval |= CleaverCUDA::kIsEvaluated;
  //There are no triples if one of its edges doesn't have a cut.
  std::array<CleaverCUDA::Edge*,3> edges =
      geos.GetFaceEdges(Idx3(i,j,k),num,(Definitions::tet_index)(-1));
  for (size_t t = 0; t < 3; t++)
    if (!(edges[t]->isCut_eval & CleaverCUDA::kIsCut)) return;
  face->isCut_eval |= CleaverCUDA::kIsCut;
  //get the 3 vertices for the face.
  //  std::array<std::array<float,3>,3> pts;
  //  for (size_t x = 0; x < 3; x++)
  //    for (size_t y = 0; y < 3; y++)
  //      pts[x][y] = edges[x]->cut_loc[y];
  //  std::array<float,3> res = FindTriplePoint(pts);
  //  for (size_t x = 0; x < 3; x++)
  //    face->cut_loc[x] = res[x];
}

void CleaverUtility::CalculateQuadruple(SimpleTet *tet,
                                        Definitions::tet_index num,
                                        size_t i, size_t j, size_t k,
                                        SimpleGeometry& geos) {
  tet->isCut_eval |= CleaverCUDA::kIsEvaluated;
  //There are no triples if one of its edges doesn't have a cut.
  std::array<SimpleFace*,4> faces =
      geos.GetTetFaces(Idx3(i,j,k),num);
  //find 3 & 4 edge cut cases.
  size_t num_cuts = 0;
  std::array<CleaverCUDA::Edge*,6> edges =
      geos.GetTetEdges(Idx3(i,j,k),num);
  for (auto a : edges) if ((a->isCut_eval & CleaverCUDA::kIsCut)) {
    num_cuts++;
    if (num_cuts > 2) {
      tet->isCut_eval |= CleaverCUDA::kHasStencil;
      break;
    }
  }
  // check for quadruple point
  for (size_t t = 0; t < 4; t++)
    if (!(faces[t]->isCut_eval & CleaverCUDA::kIsCut)) return;
  tet->isCut_eval |= CleaverCUDA::kIsCut;
}

void CleaverUtility::CalculateQuadruples(SimpleGeometry& geos) {
  clock_t start = clock();
  size_t count_quads = 0;
  for(size_t i = 0; i < w_; i ++)
    for(size_t j = 0; j < h_; j ++)
      for(size_t k = 0; k < d_; k ++) {
        if(cut_cells_[Idx3(i,j,k)]) {
          for(size_t l = 0; l < Definitions::kTetsPerCell; l++) {
            SimpleTet* tet = geos.GetTet(Idx3(i,j,k),
                                         (Definitions::tet_index)l);
            if(!(tet->isCut_eval & CleaverCUDA::kIsEvaluated)) {
              CalculateQuadruple(tet,(Definitions::tet_index)l,i,j,k,geos);
              if ((tet->isCut_eval & CleaverCUDA::kIsCut)) count_quads++;
            }
          }
        }
      }
  double duration = ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "Found all Quadruples: (" << count_quads << ")\t" <<
    duration << " sec." << std::endl;
}

void CleaverUtility::StencilFaces(SimpleGeometry& geos) {
  clock_t start = clock();
  faces_.clear(); verts_.clear();
  for(size_t i = 0; i < w_; i ++)
    for(size_t j = 0; j < h_; j ++)
      for(size_t k = 0; k < d_; k ++) {
        if(cut_cells_[Idx3(i,j,k)]) {
          for (size_t l = 0; l < Definitions::kTetsPerCell; l++) {
            SimpleTet * tet = geos.GetTet(Idx3(i,j,k),
                                          (Definitions::tet_index)l);
            //don't repeat tets
            if (!(tet->isCut_eval & CleaverCUDA::kIsStenciled)) {
              tet->isCut_eval |= CleaverCUDA::kIsStenciled;
              if ((tet->isCut_eval & CleaverCUDA::kHasStencil))
                AccumulateVertsAndFaces(
                    tet,(Definitions::tet_index)l,i,j,k,geos);
            }
          }
        }
      }
  double duration = ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "Generalized/Stenciled Faces.\t" <<
    duration << " sec." << std::endl;
}

void CleaverUtility::OutputToFile() {
  clock_t start = clock();
  std::ofstream file((output_ + ".ply"));
  file << "ply" << std::endl;
  file << "format ascii 1.0" << std::endl;
  file << "comment " << w_ << " " << h_ << " " << d_ << std::endl;
  file << "element vertex " << verts_.size() << std::endl;
  file << "property float x " << std::endl;
  file << "property float y " << std::endl;
  file << "property float z " << std::endl;
  file << "element face " << faces_.size() << std::endl;
  file << "property list uchar int vertex_index" << std::endl;
  file << "element color " << faces_.size() << std::endl;
  file << "property list uchar int face_color" << std::endl;
  file << "end_header" << std::endl;
  for(auto a : verts_)
    file << a[0] << " " << a[1] << " " << a[2] << std::endl;
  for(auto a : faces_)
    file << "3 " << a[0] << " " << a[1] << " " << a[2] << std::endl;
  for(auto a : faces_)
    file << a[3] << std::endl;
  file.close();
  double duration = ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
  if (verbose_)
    std::cout << "Wrote faces to file: " << output_ << ".ply\t" <<
    duration << " sec." << std::endl;
}

std::pair<std::vector<std::array<float,3>>,std::vector<std::array<size_t,4>>>
CleaverUtility::GetVertsFacesFromNRRD(std::vector<std::string> &files,
                                      float scales[3],
                                      std::array<size_t,3> &res,
                                      bool useScale, bool useAbs,
                                      bool useGPU, bool addAir) {
  inputs_.clear();
  for(auto a : files) inputs_.push_back(a);
  for(size_t x = 0; x < 3; x++) {
    if((scaled_resolution_ = useScale))
      scale_[x] = scales[x];
    if((absolute_resolution_ = useAbs))
      res_[x] = res[x];
  }
  this->use_GPU_ = useGPU;
  this->addAir_ = addAir;
  /***************************CPU*****************************/
  //  Load Data
  if(!LoadNRRDs()) {
    std::cerr << "Error reading data or allocating memory!" << std::endl;
  } else {
    // Find all of the max materials
    FindMaxes();
    //  Find all of the cut cells
    FindCutCells();
    //allocate memory for the geometry.
    size_t w,h,d;
    w = this->w();
    h = this->h();
    d = this->d();
    SimpleGeometry geos(w,h,d,verbose());
    if(!geos.Valid()) {
      std::cerr << "Error allocating geometric memory!" << std::endl;
    } else {
      //geos.TestCells();
      //  Calculate the cuts
      CalculateCuts(geos);
      //  Calculate the triples
      CalculateTriples(geos);
      //  Calculate the quadruples
      CalculateQuadruples(geos);
      //  Stencil Faces
      StencilFaces(geos);
    }
  }
  return std::pair<std::vector<std::array<float,3>>,
      std::vector<std::array<size_t,4>>>(verts_,faces_);
}

void CleaverUtility::GetVertsFacesFromNRRD(int argc, char *argv[]) {
  ParseInput(argc,argv);
  std::vector<std::string> ins;
  for(auto a : inputs_) ins.push_back(a);
  GetVertsFacesFromNRRD(ins,scale_,res_,
                        scaled_resolution_,absolute_resolution_,
                        use_GPU_,addAir_);
}

void CleaverUtility::AccumulateVertsAndFaces(
    SimpleTet * tet,
    Definitions::tet_index num,
    size_t i, size_t j, size_t k, SimpleGeometry& geos) {
  //set this tet to be evaluated for mesh faces.
  std::array<Definitions::tri_index,4> tet_faces_num =
      geos.GetTetFacesNum(num);
  //first case is when this tet has a quad cut (6 edge cut)
  if ((tet->isCut_eval & CleaverCUDA::kIsCut)) {
    //the first vertex is always the quadruple point
    std::array<float,3> v1 = {{0,0,0}};
    std::array<CleaverCUDA::Edge*,6> tedges =
        geos.GetTetEdges(Idx3(i,j,k),num);
    for (auto vv : tedges) {
      for (size_t x = 0; x < 3; x++)
        v1[x] += (vv->cut_loc[x]);
    }
    for (size_t x = 0; x < 3; x++) v1[x] /= 6.;
    for (size_t f = 0; f < 4; f++) {
      // the second vertex is 1 of the 4 triples points
      std::array<float,3> v2 = {{0,0,0}};
      std::array<CleaverCUDA::Edge *,3> edges =
          geos.GetFaceEdges(Idx3(i,j,k),tet_faces_num[f],num);
      for (size_t e = 0; e < 3; e++)
        for (size_t x = 0; x < 3; x++)
          v2[x] += edges[e]->cut_loc[x];
      for (size_t x = 0; x < 3; x++) v2[x] /= 3.;
      for (size_t e = 0; e < 3; e++) {
        // the third vertex is the cut point for each edge on each face.
        std::array<float,3> v3 = {{0,0,0}};
        for (size_t x = 0; x < 3; x++)
          v3[x] = edges[e]->cut_loc[x];
        AddFace({{ v1, v2, v3 }}, (edges[0]->isCut_eval >>
            CleaverCUDA::kMaterial));
      }
    }
    return;
  }
  //now look for 2 triple case (5 edge cut)
  size_t triple_count = 0;
  std::array<SimpleFace*,4> tet_faces =
      geos.GetTetFaces(Idx3(i,j,k),num);
  for (auto a : tet_faces) if ((a->isCut_eval &
      CleaverCUDA::kIsCut)) triple_count++;
  if (triple_count == 2) {
    //get the 2 triple points
    std::array<float,3> tpA, tpB, shared_cp, fa1, fa2, fb1, fb2;
    tpA = tpB = shared_cp = fa1 = fa2 = fb1 = fb2 = {{0,0,0}};
    std::array<CleaverCUDA::Edge *,3> edges1, edges2, edges3;
    edges1 = edges2 = edges3 = {{nullptr, nullptr, nullptr}};
    for (size_t x = 0; x < 3; x++) {
      edges1[x] = nullptr;
      edges2[x] = nullptr;
      edges3[x] = nullptr;
    }
    bool first = true;
    for (size_t f = 0; f < 4; f++)
      if ((tet_faces[f]->isCut_eval & CleaverCUDA::kIsCut)) {
        if (first) {
          edges1 = geos.GetFaceEdges(Idx3(i,j,k),tet_faces_num[f],num);
          for (auto aa : edges1)
            for (size_t x = 0; x < 3; x++)
              tpA[x] += aa->cut_loc[x];
          for (size_t x = 0; x < 3; x++) tpA[x] /= 3.;
        } else {
          edges2 = geos.GetFaceEdges(Idx3(i,j,k),tet_faces_num[f],num);
          for (auto aa : edges2)
            for (size_t x = 0; x < 3; x++)
              tpB[x] += aa->cut_loc[x];
          for (size_t x = 0; x < 3; x++) tpB[x] /= 3.;
        }
        if (!first && edges3[0]) break;
        first = false;
      } else if (!edges3[0]) {
        edges3 = geos.GetFaceEdges(Idx3(i,j,k),tet_faces_num[f],num);
      }
    //get the shared edge cut
    std::array<CleaverCUDA::Edge*,5> edges =
    {{nullptr,nullptr,nullptr, nullptr,nullptr}};
    for (auto& a : edges) a = nullptr;
    for (auto a: edges1) {
      for (auto b : edges2)
        if (a == b) { edges[0] = a; break; }
      if (edges[0]) break;
    }
    //get the independent edge cuts
    edges[1] = (edges[0]==edges1[0])?edges1[1]:edges1[0];
    edges[2] = (edges[0]==edges1[2])?edges1[1]:edges1[2];
    edges[3] = (edges[0]==edges2[0])?edges2[1]:edges2[0];
    edges[4] = (edges[0]==edges2[2])?edges2[1]:edges2[2];
    //determine which face each edge is on
    bool on_third[4] = { false, false, false, false };
    for (size_t t = 0; t < 4; t++)
      for (size_t u = 0; u < 3; u++)
        on_third[t] |= (edges3[u] == edges[t+1]);
    if (on_third[0] != on_third[2]) {
      CleaverCUDA::Edge * tmp = edges[1];
      edges[1] = edges[2];
      edges[2] = tmp;
    }
    for (size_t x = 0; x < 3; x++) {
      shared_cp[x] = edges[0]->cut_loc[x];
      fa1[x]       = edges[1]->cut_loc[x];
      fa2[x]       = edges[2]->cut_loc[x];
      fb1[x]       = edges[3]->cut_loc[x];
      fb2[x]       = edges[4]->cut_loc[x];
    }
    //now add the 5 faces.
    AddFace({{fa1,tpB,fb1}},
            (edges[1]->isCut_eval >> CleaverCUDA::kMaterial));
    AddFace({{fa1,tpA,tpB}},
            (edges[1]->isCut_eval >> CleaverCUDA::kMaterial));
    AddFace({{fa2,tpB,fb2}},
            (edges[2]->isCut_eval >> CleaverCUDA::kMaterial));
    AddFace({{fa2,tpA,tpB}},
            (edges[2]->isCut_eval >> CleaverCUDA::kMaterial));
    AddFace({{tpA,tpB,shared_cp}},
            (edges[0]->isCut_eval >> CleaverCUDA::kMaterial));
    return;
  }
  //now deal with 3 & 4 edge cut cases.
  std::array<CleaverCUDA::Edge *,3> edges = {{nullptr,nullptr,nullptr}};
  std::vector<std::pair<std::array<float,3>,std::vector<size_t>> > pts;
  std::vector<char> mats;
  for (auto& a : pts)
    a.first = {{0,0,0}};
  for (size_t f = 0; f < 4; f++) {
    std::array<CleaverCUDA::Edge *,3> edges =
        geos.GetFaceEdges(Idx3(i,j,k),tet_faces_num[f],num);
    for(auto a : edges)
      if ((a->isCut_eval & CleaverCUDA::kIsCut)) {
        //add the point to list of points.
        std::array<float,3> pt = {{a->cut_loc[0],
            a->cut_loc[1],a->cut_loc[2]}};
        bool found = false;
        for(size_t t = 0; t < pts.size(); t++)
          if ((pts[t].first[0] == pt[0]) &&
              (pts[t].first[1] == pt[1]) &&
              (pts[t].first[2] == pt[2])) {
            pts[t].second.push_back(f);
            found = true;
            break;
          }
        if (!found) {
          std::vector<size_t> tmp; tmp.push_back(f);
          pts.push_back(std::pair<std::array<float,3>,
                        std::vector<size_t>>(pt,tmp));
          mats.push_back(a->isCut_eval >> CleaverCUDA::kMaterial);
        }
      }
  }
  if (mats.size() == 0) return;
  //find the most common material.
  std::vector<std::pair<char,char> > mat_cnts;
  for (auto a : mats) {
    bool found = false;
    for (size_t b = 0; b < mat_cnts.size(); b++)
      if (mat_cnts[b].first == a) {
        found = true;
        mat_cnts[b].second ++;
      }
    if (!found) mat_cnts.push_back(std::pair<char,char>(a,1));
  }
  char face_mat = mat_cnts[0].first;
  char greatest = mat_cnts[0].second;
  for (size_t b = 1; b < mat_cnts.size(); b++)
    if (mat_cnts[b].second > greatest) {
      greatest = mat_cnts[b].second;
      face_mat = mat_cnts[b].first;
    }
  if (pts.size() == 4) {
    size_t C = 3;
    if ((pts[0].second[0] != pts[2].second[0]) &&
        (pts[0].second[0] != pts[2].second[1]) &&
        (pts[0].second[1] != pts[2].second[0]) &&
        (pts[0].second[1] != pts[2].second[1])) C = 2;
    size_t D = (C==2)?3:2;
    AddFace({{pts[0].first,pts[C].first,pts[1].first}},face_mat);
    AddFace({{pts[0].first,pts[C].first,pts[D].first}},face_mat);
  }
  //now deal with 3 edge cut cases.
  if (pts.size() == 3) {
    AddFace({{pts[0].first,pts[1].first,pts[2].first}},face_mat);
  }
}

void CleaverUtility::AddFace(
    std::array<std::array<float,3>,3> face, char mat) {
  //shift to center on 0,0,0.
  for (auto &a : face) {
    a[0] -= (ww_ >> 1);
    a[1] -= (hh_ >> 1);
    a[2] -= (dd_ >> 1);
  }
  std::array<size_t,4> new_face = {{0,0,0,static_cast<size_t>(mat)}};
  for (size_t i = 0; i < 3; i++) {
    new_face[i] = verts_.size();
    verts_.push_back(face[i]);
  }
  faces_.push_back(new_face);
}
