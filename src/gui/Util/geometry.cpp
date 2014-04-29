/**
 * @file geometry.cpp
 * This file contains the basic geometries and
 * geometry access for CleaverCUDA.
 * @version Feb 26, 2014
 * @author: Brig Bagley
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include "geometry.h"


SimpleGeometry::SimpleGeometry(size_t w, size_t h, size_t d, bool verbose)
: w_(w), h_(h), d_(d) {
  //allocate verts, faces, edges, and tets.
  size_t total_cells = w * h * d;
  size_t total_cells_buf = (w+1) * (h+1) * (d+1);
  //total edges (shared, no duplicates)
  num_inner_edges_ = 8 * total_cells_buf;
  inner_edges_ = new CleaverCUDA::Edge[num_inner_edges_];
  num_dual_edges_ = num_axis_edges_ = 3 * total_cells_buf;
  dual_edges_ = new CleaverCUDA::Edge[num_dual_edges_];
  axis_edges_ = new CleaverCUDA::Edge[num_axis_edges_];
  for (size_t t = 0; t < num_inner_edges_; t++)
    inner_edges_[t].isCut_eval = 0;
  for (size_t t = 0; t < num_dual_edges_; t++)
    dual_edges_[t].isCut_eval = 0;
  for (size_t t = 0; t < num_axis_edges_; t++)
    axis_edges_[t].isCut_eval = 0;
  //total faces (shared, no duplicates, inner and outer)
  num_inner_faces_ = num_outer_faces_ = num_tets_ = 12 * total_cells_buf;
  inner_faces_ = new SimpleFace[num_inner_faces_];
  outer_faces_ = new SimpleFace[num_outer_faces_];
  tets_ = new SimpleTet[num_tets_];
  size_t memory_footprint =
      (sizeof(float) * total_cells +
          sizeof(bool) * total_cells +
          sizeof(char) * total_cells +
          num_inner_edges_ * sizeof(CleaverCUDA::Edge) +
          num_dual_edges_ * sizeof(CleaverCUDA::Edge) +
          num_axis_edges_ * sizeof(CleaverCUDA::Edge) +
          num_inner_faces_ * sizeof(SimpleFace) +
          num_outer_faces_ * sizeof(SimpleFace) +
          num_tets_ * sizeof(SimpleTet)) >> 20;
  if (!Valid())
    std::cerr << "Failed to allocate memory for geometry." << std::endl;
  if(verbose)
    std::cout << "Memory:\t\t\t\t" << memory_footprint << " MB (" <<
    w_ << " x " << h_ << " x " << d_ << ")" << std::endl;
}

SimpleGeometry::~SimpleGeometry() {
  delete[] inner_edges_;
  delete[] dual_edges_;
  delete[] axis_edges_;
  delete[] inner_faces_;
  delete[] outer_faces_;
  delete[] tets_;
}

bool SimpleGeometry::Valid() {
  return inner_edges_ && dual_edges_ &&
      axis_edges_ && inner_faces_ && outer_faces_ && tets_;
}

std::array<std::array<float,3>,3> SimpleGeometry::GetFaceVertices(
    Definitions::tri_index num,
    size_t i, size_t j, size_t k,
    std::array<float,3> scale) {
  float x = static_cast<float>(i) * scale[0];
  float y = static_cast<float>(j) * scale[1];
  float z = static_cast<float>(k) * scale[2];
  std::array<float,3> v1 = {{x,y,z}};
  std::array<float,3> v2 = {{x+scale[0],y,z}};
  std::array<float,3> v3 = {{x,y+scale[1],z}};
  std::array<float,3> v4 = {{x+scale[0],y+scale[1],z}};
  std::array<float,3> v5 = {{x,y,z+scale[2]}};
  std::array<float,3> v6 = {{x+scale[0],y,z+scale[2]}};
  std::array<float,3> v7 = {{x,y+scale[1],z+scale[2]}};
  std::array<float,3> v8 = {{x+scale[0],y+scale[1],z+scale[2]}};
  std::array<float,3> v9 = {{x+.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v10 = {{x-.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v11 = {{x+1.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v12 = {{x+.5f * scale[0],
      y-.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v13 = {{x+.5f * scale[0],
      y+1.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v14 = {{x+.5f * scale[0],
      y+.5f * scale[1],z-.5f * scale[2]}};
  std::array<float,3> v15 = {{x+.5f * scale[0],
      y+.5f * scale[1],z+1.5f * scale[2]}};
  std::array<std::array<float,3>,3> ans;
  switch (num) {
    //inner faces
    case Definitions::FUL: ans[0] = v3; ans[1] = v9; ans[2] = v7; break;
    case Definitions::FUR: ans[0] = v4; ans[1] = v9; ans[2] = v8; break;
    case Definitions::FUF: ans[0] = v3; ans[1] = v9; ans[2] = v4; break;
    case Definitions::FUB: ans[0] = v7; ans[1] = v9; ans[2] = v8; break;
    case Definitions::FLL: ans[0] = v1; ans[1] = v9; ans[2] = v5; break;
    case Definitions::FLR: ans[0] = v2; ans[1] = v9; ans[2] = v6; break;
    case Definitions::FLF: ans[0] = v1; ans[1] = v9; ans[2] = v2; break;
    case Definitions::FLB: ans[0] = v5; ans[1] = v9; ans[2] = v6; break;
    case Definitions::FFL: ans[0] = v1; ans[1] = v9; ans[2] = v3; break;
    case Definitions::FFR: ans[0] = v2; ans[1] = v9; ans[2] = v4; break;
    case Definitions::FBL: ans[0] = v5; ans[1] = v9; ans[2] = v7; break;
    case Definitions::FBR: ans[0] = v6; ans[1] = v9; ans[2] = v8; break;
    //left faces
    case Definitions::FLUF: ans[0] = v10; ans[1] = v3; ans[2] = v9; break;
    case Definitions::FLUB: ans[0] = v10; ans[1] = v7; ans[2] = v9; break;
    case Definitions::FLLF: ans[0] = v10; ans[1] = v1; ans[2] = v9; break;
    case Definitions::FLLB: ans[0] = v10; ans[1] = v5; ans[2] = v9; break;
    //right faces
    case Definitions::FRUF: ans[0] = v9; ans[1] = v4; ans[2] = v11; break;
    case Definitions::FRUB: ans[0] = v9; ans[1] = v8; ans[2] = v11; break;
    case Definitions::FRLF: ans[0] = v9; ans[1] = v2; ans[2] = v11; break;
    case Definitions::FRLB: ans[0] = v9; ans[1] = v6; ans[2] = v11; break;
    //up faces
    case Definitions::FUFL: ans[0] = v9; ans[1] = v3; ans[2] = v13; break;
    case Definitions::FUBL: ans[0] = v9; ans[1] = v7; ans[2] = v13; break;
    case Definitions::FUFR: ans[0] = v9; ans[1] = v4; ans[2] = v13; break;
    case Definitions::FUBR: ans[0] = v9; ans[1] = v8; ans[2] = v13; break;
    //down faces
    case Definitions::FDFL: ans[0] = v12; ans[1] = v1; ans[2] = v9; break;
    case Definitions::FDBL: ans[0] = v12; ans[1] = v5; ans[2] = v9; break;
    case Definitions::FDFR: ans[0] = v12; ans[1] = v2; ans[2] = v9; break;
    case Definitions::FDBR: ans[0] = v12; ans[1] = v6; ans[2] = v9; break;
    //front faces
    case Definitions::FFUL: ans[0] = v14; ans[1] = v3; ans[2] = v9; break;
    case Definitions::FFUR: ans[0] = v14; ans[1] = v4; ans[2] = v9; break;
    case Definitions::FFLL: ans[0] = v14; ans[1] = v1; ans[2] = v9; break;
    case Definitions::FFLR: ans[0] = v14; ans[1] = v2; ans[2] = v9; break;
    //back faces
    case Definitions::FBUL: ans[0] = v9; ans[1] = v7; ans[2] = v15; break;
    case Definitions::FBUR: ans[0] = v9; ans[1] = v8; ans[2] = v15; break;
    case Definitions::FBLL: ans[0] = v9; ans[1] = v5; ans[2] = v15; break;
    case Definitions::FBLR: ans[0] = v9; ans[1] = v6; ans[2] = v15; break;
  }
  return ans;
}

std::array<std::array<float,3>,4> SimpleGeometry::GetTetVertices(
    Definitions::tet_index num,
    size_t i, size_t j, size_t k,
    std::array<float,3> scale) {
  float x = static_cast<float>(i) * scale[0];
  float y = static_cast<float>(j) * scale[1];
  float z = static_cast<float>(k) * scale[2];
  std::array<float,3> v1 = {{x,y,z}};
  std::array<float,3> v2 = {{x+scale[0],y,z}};
  std::array<float,3> v3 = {{x,y+scale[1],z}};
  std::array<float,3> v4 = {{x+scale[0],y+scale[1],z}};
  std::array<float,3> v5 = {{x,y,z+scale[2]}};
  std::array<float,3> v6 = {{x+scale[0],y,z+scale[2]}};
  std::array<float,3> v7 = {{x,y+scale[1],z+scale[2]}};
  std::array<float,3> v8 = {{x+scale[0],y+scale[1],z+scale[2]}};
  std::array<float,3> v9 = {{x+.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v10 = {{x-.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v11 = {{x+1.5f * scale[0],
      y+.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v12 = {{x+.5f * scale[0],
      y-.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v13 = {{x+.5f * scale[0],
      y+1.5f * scale[1],z+.5f * scale[2]}};
  std::array<float,3> v14 = {{x+.5f * scale[0],
      y+.5f * scale[1],z-.5f * scale[2]}};
  std::array<float,3> v15 = {{x+.5f * scale[0],
      y+.5f * scale[1],z+1.5f * scale[2]}};
  std::array<std::array<float,3>,4> ans;
  switch (num) {
    //left tets
    case Definitions::TLU:
      ans[0] = v3; ans[1] = v7; ans[2] = v10; ans[3] = v9; break;
    case Definitions::TLL:
      ans[0] = v1; ans[1] = v5; ans[2] = v10; ans[3] = v9; break;
    case Definitions::TLF:
      ans[0] = v1; ans[1] = v3; ans[2] = v10; ans[3] = v9; break;
    case Definitions::TLB:
      ans[0] = v5; ans[1] = v7; ans[2] = v10; ans[3] = v9; break;
      //right tets
    case Definitions::TRU:
      ans[0] = v4; ans[1] = v8; ans[2] = v9; ans[3] = v11; break;
    case Definitions::TRL:
      ans[0] = v2; ans[1] = v6; ans[2] = v9; ans[3] = v11; break;
    case Definitions::TRF:
      ans[0] = v2; ans[1] = v4; ans[2] = v9; ans[3] = v11; break;
    case Definitions::TRB:
      ans[0] = v6; ans[1] = v8; ans[2] = v9; ans[3] = v11; break;
      //front tets
    case Definitions::TFT:
      ans[0] = v3; ans[1] = v4; ans[2] = v14; ans[3] = v9; break;
    case Definitions::TFB:
      ans[0] = v1; ans[1] = v2; ans[2] = v14; ans[3] = v9; break;
    case Definitions::TFL:
      ans[0] = v1; ans[1] = v3; ans[2] = v14; ans[3] = v9; break;
    case Definitions::TFR:
      ans[0] = v2; ans[1] = v4; ans[2] = v14; ans[3] = v9; break;
      //back tets
    case Definitions::TBT:
      ans[0] = v7; ans[1] = v8; ans[2] = v9; ans[3] = v15; break;
    case Definitions::TBB:
      ans[0] = v5; ans[1] = v6; ans[2] = v9; ans[3] = v15; break;
    case Definitions::TBL:
      ans[0] = v5; ans[1] = v7; ans[2] = v9; ans[3] = v15; break;
    case Definitions::TBR:
      ans[0] = v6; ans[1] = v8; ans[2] = v9; ans[3] = v15; break;
      //down tets
    case Definitions::TDF:
      ans[0] = v1; ans[1] = v2; ans[2] = v12; ans[3] = v9; break;
    case Definitions::TDB:
      ans[0] = v5; ans[1] = v6; ans[2] = v12; ans[3] = v9; break;
    case Definitions::TDL:
      ans[0] = v1; ans[1] = v5; ans[2] = v12; ans[3] = v9; break;
    case Definitions::TDR:
      ans[0] = v2; ans[1] = v6; ans[2] = v12; ans[3] = v9; break;
      //down tets
    case Definitions::TUF:
      ans[0] = v3; ans[1] = v4; ans[2] = v9; ans[3] = v13; break;
    case Definitions::TUB:
      ans[0] = v7; ans[1] = v8; ans[2] = v9; ans[3] = v13; break;
    case Definitions::TUL:
      ans[0] = v3; ans[1] = v7; ans[2] = v9; ans[3] = v13; break;
    case Definitions::TUR:
      ans[0] = v4; ans[1] = v8; ans[2] = v9; ans[3] = v13; break;
  }
  return ans;
}

std::array<CleaverCUDA::Edge*,3> SimpleGeometry::GetFaceEdges(
    size_t cell, Definitions::tri_index fnum, Definitions::tet_index tnum) {
  std::array<CleaverCUDA::edge_index,3> edges = GetFaceEdgesNum(fnum);
  int64_t first = cell, second = cell, cel = cell;
  if (tnum < 0 || tnum >23) {
    if (fnum < 12)
      first = second = cel;
    else if (fnum < 16)
      first = cel - 1;
    else if (fnum < 20)
      second = cel + 1;
    else if (fnum < 24)
      second = cel + w_;
    else if (fnum < 28)
      first = cel - w_;
    else if (fnum < 32)
      first = cel - w_*h_;
    else
      second = cel + w_*h_;
  } else if (fnum > 11) {
    if (tnum < 4 && (11 < fnum && fnum < 16) )
      first = cel - 1;
    else if (3 < tnum && tnum < 8 && (15 < fnum && fnum < 20) )
      second = cel + 1;
    else if (7 < tnum && tnum < 12 && (27 < fnum && fnum < 32) )
      first = cel - w_*h_;
    else if (11 < tnum && tnum < 16 && (31 < fnum && fnum < 36) )
      second = cel + w_*h_;
    else if (15 < tnum && tnum < 20 && (23 < fnum && fnum < 28))
      first = cel - w_;
    else if (19 < tnum && tnum < 24 &&  (19 < fnum && fnum < 24) )
      second = cel + w_;
  } else {
    if (tnum < 4 &&
        (fnum == 1 || fnum == 5 || fnum == 9 || fnum == 11))
      first = second = cel = cel - 1;
    else if (3 < tnum && tnum < 8 &&
        (fnum == 0 || fnum == 4 || fnum == 8 || fnum == 10))
      first = second = cel = cel + 1;
    else if (7 < tnum && tnum < 12 &&
        (fnum == 3 || fnum == 7 || fnum == 10 || fnum == 11))
      first = second = cel = cel - w_*h_;
    else if (11 < tnum && tnum < 16 &&
        (fnum == 2 || fnum == 6 || fnum == 8 || fnum == 9))
      first = second = cel = cel + w_*h_;
    else if (15 < tnum && tnum < 20 && fnum < 4)
      first = second = cel = cel - w_;
    else if (19 < tnum && tnum < 24 && 3 < fnum && fnum < 8)
      first = second = cel = cel + w_;
  }
  std::array<CleaverCUDA::Edge*,3> out = {{
      GetEdge(first,edges[0]),
      GetEdge(second,edges[1]),
      GetEdge(cel,edges[2]) }};
  return out;
}

std::array<CleaverCUDA::edge_index,3> SimpleGeometry::GetFaceEdgesNum(
    Definitions::tri_index n) {
  std::array<size_t,3> edges = {{0,0,0}};
  size_t num = static_cast<size_t>(n);
  bool even   = ((num % 2) == 0);
  bool gto    = ((num % 4) > 1);
  size_t add  = ((num % 4) == 1)?2:1;
  size_t add2 = gto?2:1;
  size_t add4 = gto?4:0;
  if (num < 4)
    edges = {{ (even?0:add), add2+(even?0:add), num+14 }};
  else if (num < 8)
    edges = {{ 4+(even?0:add), 4+add2+(even?0:add), num+14 }};
  else if (num < 12)
    edges = {{ 4lu + (even?0lu:2lu) + (gto?1lu:0lu),
        (even?0lu:2lu) + (gto?1lu:0lu), num+14lu }};
  else if (num < 20)
    edges = {{ (2lu+add4+(even?0lu:1lu)),
        (add4+(even?0lu:1lu)), ((num<16)?8lu:9lu) }};
  else if (num < 28)
    edges = {{ (num%4lu), ((num%4lu) + 4lu), ((num<24lu)?10lu:11lu) }};
  else
    edges = {{ (((num%4)*2lu)+1lu),
        ((num%4lu)*2lu), ((num<32lu)?12lu:13lu) }};
  std::array<CleaverCUDA::edge_index,3> e;
  for (size_t a = 0; a < 3; a++)
    e[a] = static_cast<CleaverCUDA::edge_index>(edges[a]);
  return e;
}

std::array<SimpleFace*,4> SimpleGeometry::GetTetFaces(
    size_t cell, Definitions::tet_index n) {
  std::array<Definitions::tri_index,4> faces = GetTetFacesNum(n);
  int64_t first = cell, second = cell, cel = cell;
  if (n < 4)
    first = cel - 1;
  else if (n < 8)
    second = cel + 1;
  else if (n < 12)
    first = cel - w_ * h_;
  else if (n < 16)
    second = cel + w_ * h_;
  else if (n < 20)
    first = cel - w_;
  else
    second = cel + w_;
  std::array<SimpleFace*,4> out = {{ GetFace(first,faces[0]),
      GetFace(second,faces[1]), GetFace(cel,faces[2]),
      GetFace(cel,faces[3]) }};
  return out;
}

std::array<CleaverCUDA::Edge*,6> SimpleGeometry::GetTetEdges(
    size_t cell, Definitions::tet_index n) {
  std::array<CleaverCUDA::Edge*,6> ans;
  std::array<Definitions::tri_index,4> fcs = GetTetFacesNum(n);
  size_t idx = 0;
  for (auto a : fcs) {
    std::array<CleaverCUDA::Edge*,3> edges = GetFaceEdges(cell,a,n);
    for (auto b: edges) {
      bool found = false;
      for (size_t x = 0; x < idx; x++)
        found |= (ans[x] == b);
      if (!found) {
        ans[idx] = b;
        idx++;
        if(idx >= 6) return ans;
      }
    }
  }
  return ans;
}

std::array<Definitions::tri_index,4> SimpleGeometry::GetTetFacesNum(
    Definitions::tet_index n) {
  std::array<size_t,4> faces = {{0,0,0,0}};
  size_t num = static_cast<size_t>(n);

  bool even   = ((num % 2) == 0);
  bool fh    = ((num % 4) < 2);
  bool fs    = ((num % 8) < 4);

  if (num < 8)
    faces = {{ fh?(even?1ul:5ul):(even?9ul:11ul),
        fh?(even?0ul:4ul):(even?8ul:10ul),
            fh?(even?(fs?12ul:16ul):(fs?14ul:18ul)):
                (even?(fs?14ul:18ul):(fs?15ul:19ul)),
                fh?(even?(fs?13ul:17ul):(fs?15ul:19ul)):
                    (even?(fs?12ul:16ul):(fs?13ul:17ul)) }};
  else if (num < 16)
    faces = {{ fh?(even?3ul:7ul):(even?10ul:11ul),
        fh?(even?2ul:6ul):(even?8ul:9ul),
            fh?(even?(fs?28ul:32ul):(fs?30ul:34ul)):
                (even?(fs?30ul:34ul):(fs?31ul:35ul)),
                fh?(even?(fs?29ul:33ul):(fs?31ul:35ul)):
                    (even?(fs?28ul:32ul):(fs?29ul:33ul)) }};
  else
    faces = {{ fh?(even?2ul:3ul):(even?0ul:1ul),
        fh?(even?6ul:7ul):(even?4ul:5ul),
            fh?(even?(fs?24ul:20ul):(fs?25ul:21ul)):
                (even?(fs?24ul:20ul):(fs?26ul:22ul)),
                fh?(even?(fs?26ul:22ul):(fs?27ul:23ul)):
                    (even?(fs?25ul:21ul):(fs?27ul:23ul)) }};
  std::array<Definitions::tri_index,4> f;
  for (size_t a = 0; a < 4; a++)
    f[a] = static_cast<Definitions::tri_index>(faces[a]);
  return f;
}

size_t SimpleGeometry::GetInnerFaceIdx(size_t cell, size_t face) {
  return cell * 12 + face;
}

size_t SimpleGeometry::GetOuterFaceIdx(size_t cell,
                                       Definitions::cell_index dir,
                                       unsigned char num){
  return OffsetAlgorithm(cell,dir,4,num);
}

size_t SimpleGeometry::GetInnerEdgeIdx(size_t cell, size_t edge) {
  return cell * 8 + edge;
}

size_t SimpleGeometry::GetDualEdgeIdx(size_t cell,
                                      CleaverCUDA::edge_index edge){
  size_t offset = w_*h_*3, num = 0;
  switch (edge) {
    case CleaverCUDA::CL: offset -= 3;       num = 0; break;
    case CleaverCUDA::CR: offset -= 0;       num = 0; break;
    case CleaverCUDA::CU: offset -= 0;       num = 1; break;
    case CleaverCUDA::CD: offset -= w_*3;    num = 1; break;
    case CleaverCUDA::CF: offset -= w_*h_*3; num = 2; break;
    case CleaverCUDA::CB: offset -= 0;       num = 2; break;
    default:              offset -= 0;       num = 0; break;
  }
  return cell * 3 + offset + num;
}

size_t SimpleGeometry::GetAxisEdgeIdx(size_t cell,
                                      CleaverCUDA::edge_index edge){

  size_t offset = 0, num = 0;
  switch (edge) {
    case CleaverCUDA::UL: offset = w_;      num = 2; break;
    case CleaverCUDA::UR: offset = w_+1;    num = 2; break;
    case CleaverCUDA::UF: offset = w_;      num = 0; break;
    case CleaverCUDA::UB: offset = w_*h_+w_;num = 0; break;
    case CleaverCUDA::LL: offset = 0;       num = 2; break;
    case CleaverCUDA::LR: offset = 1;       num = 2; break;
    case CleaverCUDA::LF: offset = 0;       num = 0; break;
    case CleaverCUDA::LB: offset = w_*h_;   num = 0; break;
    case CleaverCUDA::FL: offset = 0;       num = 1; break;
    case CleaverCUDA::FR: offset = 1;       num = 1; break;
    case CleaverCUDA::BL: offset = w_*h_;   num = 1; break;
    case CleaverCUDA::BR: offset = w_*h_+1; num = 1; break;
    default: offset = w_; num = 0; break;
  }
  return cell * 3 + offset*3 + num;
}

size_t SimpleGeometry::OffsetAlgorithm(size_t cell,
                                       Definitions::cell_index dir,
                                       unsigned char multiple,
                                       unsigned char num) {
  int offset = 0;
  unsigned char multiple_t3 = multiple * 3;
  unsigned char multiple_t2 = multiple * 2;
  switch (dir) {
    case Definitions::RIGHT:
      offset = 0;
      break;
    case Definitions::LEFT:
      offset = -multiple_t3;
      break;
    case Definitions::UP:
      offset = multiple;
      break;
    case Definitions::DOWN:
      offset = -w_*multiple_t3 + multiple;
      break;
    case Definitions::BACK:
      offset = multiple_t2;
      break;
    case Definitions::FRONT:
      offset = -w_*h_*multiple_t3 + multiple_t2;
      break;
    default: break;
  }
  return cell*multiple_t3 + w_*h_*multiple_t3 - multiple_t2 + offset + num;
}

CleaverCUDA::Edge * SimpleGeometry::GetEdge(int64_t cell,
                                     CleaverCUDA::edge_index edge) {
  if (edge < 8) {
    cell += w_ * h_;
    return &inner_edges_[GetInnerEdgeIdx(cell,(size_t)edge)];
  } else if (edge < 14) {
    cell += w_ * h_;
    return &dual_edges_[GetDualEdgeIdx(cell,edge)];
  } else {
    return &axis_edges_[GetAxisEdgeIdx(cell,edge)];
  }
}

std::array<CleaverCUDA::Edge *,3> SimpleGeometry::GetEdgePointers() {
  return {{inner_edges_,dual_edges_,axis_edges_}};
}

std::array<size_t,3> SimpleGeometry::GetEdgePointersSize() {
  return {{num_inner_edges_,num_dual_edges_,num_axis_edges_}};
}

SimpleFace * SimpleGeometry::GetFace(int64_t cell,
                                     Definitions::tri_index face) {
  cell += w_ * h_;
  if (face < 12) {
    return &inner_faces_[GetInnerFaceIdx(cell,(size_t)face)];
  }
  Definitions::cell_index dir = Definitions::LEFT;
  unsigned char num = 0;
  switch (face) {
    //left face cutting faces
    case Definitions::FLUF: dir = Definitions::LEFT; num = 1; break;
    case Definitions::FLUB: dir = Definitions::LEFT; num = 0; break;
    case Definitions::FLLF: dir = Definitions::LEFT; num = 2; break;
    case Definitions::FLLB: dir = Definitions::LEFT; num = 3; break;
    //right face cutting faces
    case Definitions::FRUF: dir = Definitions::RIGHT; num = 1; break;
    case Definitions::FRUB: dir = Definitions::RIGHT; num = 0; break;
    case Definitions::FRLF: dir = Definitions::RIGHT; num = 2; break;
    case Definitions::FRLB: dir = Definitions::RIGHT; num = 3; break;
    //front face cutting faces
    case Definitions::FFUL: dir = Definitions::FRONT; num = 0; break;
    case Definitions::FFUR: dir = Definitions::FRONT; num = 1; break;
    case Definitions::FFLL: dir = Definitions::FRONT; num = 3; break;
    case Definitions::FFLR: dir = Definitions::FRONT; num = 2; break;
    //back face cutting faces
    case Definitions::FBUL: dir = Definitions::BACK; num = 0; break;
    case Definitions::FBUR: dir = Definitions::BACK; num = 1; break;
    case Definitions::FBLL: dir = Definitions::BACK; num = 3; break;
    case Definitions::FBLR: dir = Definitions::BACK; num = 2; break;
    //upper face cutting faces
    case Definitions::FUFL: dir = Definitions::UP; num = 3; break;
    case Definitions::FUFR: dir = Definitions::UP; num = 2; break;
    case Definitions::FUBL: dir = Definitions::UP; num = 0; break;
    case Definitions::FUBR: dir = Definitions::UP; num = 1; break;
    //lower face cutting faces
    case Definitions::FDFL: dir = Definitions::DOWN; num = 3; break;
    case Definitions::FDFR: dir = Definitions::DOWN; num = 2; break;
    case Definitions::FDBL: dir = Definitions::DOWN; num = 0; break;
    case Definitions::FDBR: dir = Definitions::DOWN; num = 1; break;
    default: break;
  }
  return &outer_faces_[GetOuterFaceIdx(cell,dir,num)];
}

SimpleTet * SimpleGeometry::GetTet(int64_t cell,
                                   Definitions::tet_index tet) {
  cell += w_ * h_;
  Definitions::cell_index dir = Definitions::LEFT;
  unsigned char num = 0;
  switch (tet) {
    //left face cutting faces
    case Definitions::TLU: dir = Definitions::LEFT; num = 0; break;
    case Definitions::TLL: dir = Definitions::LEFT; num = 2; break;
    case Definitions::TLF: dir = Definitions::LEFT; num = 1; break;
    case Definitions::TLB: dir = Definitions::LEFT; num = 3; break;
    //right face cutting faces
    case Definitions::TRU: dir = Definitions::RIGHT; num = 0; break;
    case Definitions::TRL: dir = Definitions::RIGHT; num = 2; break;
    case Definitions::TRF: dir = Definitions::RIGHT; num = 1; break;
    case Definitions::TRB: dir = Definitions::RIGHT; num = 3; break;
    //front face cutting faces
    case Definitions::TFT: dir = Definitions::FRONT; num = 0; break;
    case Definitions::TFB: dir = Definitions::FRONT; num = 2; break;
    case Definitions::TFL: dir = Definitions::FRONT; num = 3; break;
    case Definitions::TFR: dir = Definitions::FRONT; num = 1; break;
    //back face cutting faces
    case Definitions::TBT: dir = Definitions::BACK; num = 0; break;
    case Definitions::TBB: dir = Definitions::BACK; num = 2; break;
    case Definitions::TBL: dir = Definitions::BACK; num = 3; break;
    case Definitions::TBR: dir = Definitions::BACK; num = 1; break;
    //lower face cutting faces
    case Definitions::TDF: dir = Definitions::DOWN; num = 2; break;
    case Definitions::TDB: dir = Definitions::DOWN; num = 0; break;
    case Definitions::TDL: dir = Definitions::DOWN; num = 3; break;
    case Definitions::TDR: dir = Definitions::DOWN; num = 1; break;
    //upper face cutting faces
    case Definitions::TUF: dir = Definitions::UP; num = 2; break;
    case Definitions::TUB: dir = Definitions::UP; num = 0; break;
    case Definitions::TUL: dir = Definitions::UP; num = 3; break;
    case Definitions::TUR: dir = Definitions::UP; num = 1; break;
    default: break;
  }
  return &tets_[OffsetAlgorithm(cell,dir,4,num)];
}

std::array<size_t,3> SimpleGeometry::CreateArray(
    size_t i, size_t j, size_t k) {
  return std::array<size_t,3>({{i,j,k}});
}

void SimpleGeometry::TestCells() {
  try {
    size_t total_cells = w_*h_*d_;
    size_t w = w_, h = h_;
    // FOR ALL CELLS
    for(size_t i = 0; i < total_cells; i++) {
      ///////////////TESTING EDGES
      //inner diagonal edges must exist
      for (size_t j = 0; j < 8; j++)
        if (!GetEdge(i,(CleaverCUDA::edge_index)j))
          std::cerr << "Edge problem at cell: " << i << "Edge: " << j << std::endl;
      //Test ALL positive dual edges
      if (GetEdge(i,CleaverCUDA::CR) !=
          GetEdge(i+1,CleaverCUDA::CL))
        std::cerr << "CR/CL Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::CU) !=
          GetEdge(i+w,CleaverCUDA::CD))
        std::cerr << "CU/CD Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::CB) !=
          GetEdge(i+w*h,CleaverCUDA::CF))
        std::cerr << "CB/CF Edge problem at : " << i << std::endl;
      // Test all negative Dual edges if possible
      if (i >= 1 && GetEdge(i,CleaverCUDA::CL) !=
          GetEdge(i-1,CleaverCUDA::CR))
        std::cerr << "CL/CR Edge problem at : " << i << std::endl;
      if (i >= w && GetEdge(i,CleaverCUDA::CD) !=
          GetEdge(i-w,CleaverCUDA::CU))
        std::cerr << "CR/CL Edge problem at : " << i << std::endl;
      if (i >= w*h && GetEdge(i,CleaverCUDA::CF) !=
          GetEdge(i-w*h,CleaverCUDA::CB))
        std::cerr << "CB/CF Edge problem at : " << i << std::endl;
      //TEST ALL AXIS EDGES
      // upper left
      if (i >= 1 && GetEdge(i,CleaverCUDA::UL) !=
          GetEdge(i-1,CleaverCUDA::UR))
        std::cerr << "UL/UR Edge problem at : " << i << std::endl;
      if (i >= 1-w && GetEdge(i,CleaverCUDA::UL) !=
          GetEdge(i-1+w,CleaverCUDA::LR))
        std::cerr << "UL/LR Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::UL) !=
          GetEdge(i+w,CleaverCUDA::LL))
        std::cerr << "UL/LL Edge problem at : " << i << std::endl;
      // upper right
      if (GetEdge(i,CleaverCUDA::UR) !=
          GetEdge(i+1,CleaverCUDA::UL))
        std::cerr << "UR/UL Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::UR) !=
          GetEdge(i+w,CleaverCUDA::LR))
        std::cerr << "UR/LR Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::UR) !=
          GetEdge(i+w+1,CleaverCUDA::LL))
        std::cerr << "UR/LL Edge problem at : " << i << std::endl;
      //upper front
      if (GetEdge(i,CleaverCUDA::UF) !=
          GetEdge(i+w,CleaverCUDA::LF))
        std::cerr << "UF/LF Edge problem at : " << i << std::endl;
      if (i >= w*h && GetEdge(i,CleaverCUDA::UF) !=
          GetEdge(i-w*h,CleaverCUDA::UB))
        std::cerr << "UF/UB Edge problem at : " << i << std::endl;
      if (i>=w*h-w && GetEdge(i,CleaverCUDA::UF) !=
          GetEdge(i+w-w*h,CleaverCUDA::LB))
        std::cerr << "UF/LF Edge problem at : " << i << std::endl;
      //upper back
      if (GetEdge(i,CleaverCUDA::UB) !=
          GetEdge(i+w*h,CleaverCUDA::UF))
        std::cerr << "UB/UF Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::UB) !=
          GetEdge(i+w*h+w,CleaverCUDA::LF))
        std::cerr << "UB/LF Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::UB) !=
          GetEdge(i+w,CleaverCUDA::LB))
        std::cerr << "UB/LB Edge problem at : " << i << std::endl;
      //lower left
      if (i>=1 && GetEdge(i,CleaverCUDA::LL) !=
          GetEdge(i-1,CleaverCUDA::LR))
        std::cerr << "LL/LR Edge problem at : " << i << std::endl;
      if (i>=1+w && GetEdge(i,CleaverCUDA::LL) !=
          GetEdge(i-1-w,CleaverCUDA::UR))
        std::cerr << "LL/UR Edge problem at : " << i << std::endl;
      if (i>=w && GetEdge(i,CleaverCUDA::LL) !=
          GetEdge(i-w,CleaverCUDA::UL))
        std::cerr << "LL/UL Edge problem at : " << i << std::endl;
      //lower right
      if (GetEdge(i,CleaverCUDA::LR) !=
          GetEdge(i+1,CleaverCUDA::LL))
        std::cerr << "LR/LL Edge problem at : " << i << std::endl;
      if (i>=w && GetEdge(i,CleaverCUDA::LR) !=
          GetEdge(i-w,CleaverCUDA::UR))
        std::cerr << "LR/UR Edge problem at : " << i << std::endl;
      if (i>=w-1 && GetEdge(i,CleaverCUDA::LR) !=
          GetEdge(i-w+1,CleaverCUDA::UL))
        std::cerr << "LR/UL Edge problem at : " << i << std::endl;
      //lower front
      if (i>=w*h && GetEdge(i,CleaverCUDA::LF) !=
          GetEdge(i-w*h,CleaverCUDA::LB))
        std::cerr << "LF/LB Edge problem at : " << i << std::endl;
      if (i>=w*h+w && GetEdge(i,CleaverCUDA::LF) !=
          GetEdge(i-w*h-w,CleaverCUDA::UB))
        std::cerr << "LF/UB Edge problem at : " << i << std::endl;
      if (i>=w && GetEdge(i,CleaverCUDA::LF) !=
          GetEdge(i-w,CleaverCUDA::UF))
        std::cerr << "LF/LB Edge problem at : " << i << std::endl;
      //lower back
      if (i>=w && GetEdge(i,CleaverCUDA::LB) !=
          GetEdge(i-w,CleaverCUDA::UB))
        std::cerr << "LB/UB Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::LB) !=
          GetEdge(i+w*h,CleaverCUDA::LF))
        std::cerr << "LB/LF Edge problem at : " << i << std::endl;
      if (i>=w-w*h && GetEdge(i,CleaverCUDA::LB) !=
          GetEdge(i-w+w*h,CleaverCUDA::UF))
        std::cerr << "LB/UF Edge problem at : " << i << std::endl;
      //front left
      if (i>=1 && GetEdge(i,CleaverCUDA::FL) !=
          GetEdge(i-1,CleaverCUDA::FR))
        std::cerr << "FL/FR Edge problem at : " << i << std::endl;
      if (i>=1+w*h && GetEdge(i,CleaverCUDA::FL) !=
          GetEdge(i-1-w*h,CleaverCUDA::BR))
        std::cerr << "FL/BR Edge problem at : " << i << std::endl;
      if (i>=w*h && GetEdge(i,CleaverCUDA::FL) !=
          GetEdge(i-w*h,CleaverCUDA::BL))
        std::cerr << "FL/BL Edge problem at : " << i << std::endl;
      //front right
      if (GetEdge(i,CleaverCUDA::FR) !=
          GetEdge(i+1,CleaverCUDA::FL))
        std::cerr << "FR/FL Edge problem at : " << i << std::endl;
      if (i>=w*h && GetEdge(i,CleaverCUDA::FR) !=
          GetEdge(i-w*h,CleaverCUDA::BR))
        std::cerr << "FR/BR Edge problem at : " << i << std::endl;
      if (i>=w*h-1 && GetEdge(i,CleaverCUDA::FR) !=
          GetEdge(i-w*h+1,CleaverCUDA::BL))
        std::cerr << "FR/BL Edge problem at : " << i << std::endl;
      //back left
      if (i>=1 && GetEdge(i,CleaverCUDA::BL) !=
          GetEdge(i-1,CleaverCUDA::BR))
        std::cerr << "BL/BR Edge problem at : " << i << std::endl;
      if (i>=1-w*h && GetEdge(i,CleaverCUDA::BL) !=
          GetEdge(i-1+w*h,CleaverCUDA::FR))
        std::cerr << "BL/FR Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::BL) !=
          GetEdge(i+w*h,CleaverCUDA::FL))
        std::cerr << "BL/FL Edge problem at : " << i << std::endl;
      //back right
      if (GetEdge(i,CleaverCUDA::BR) !=
          GetEdge(i+1,CleaverCUDA::BL))
        std::cerr << "BR/BL Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::BR) !=
          GetEdge(i+1+w*h,CleaverCUDA::FL))
        std::cerr << "BR/FL Edge problem at : " << i << std::endl;
      if (GetEdge(i,CleaverCUDA::BR) !=
          GetEdge(i+w*h,CleaverCUDA::FR))
        std::cerr << "BR/FR Edge problem at : " << i << std::endl;
      //////////////TESTING FACES
      for (size_t j = 0; j < 4; j++) {
        //test inner faces
        for (size_t k = 0; k < 3; k++)
          if(!GetFace(i,(Definitions::tri_index)(j + k*4)))
            std::cerr << "Inner Face problem at : " << i <<
            ", " << (j + k*4) << std::endl;
        //test outer left
        if(i-1 >= 0 &&
            GetFace(i,(Definitions::tri_index)(j+12)) !=
                GetFace(i-1,(Definitions::tri_index)(j+16)))
          std::cerr << "Outer L/R problem at : " << i <<
          ", " << j << std::endl;
        //test outer right
        if(GetFace(i,(Definitions::tri_index)(j+16)) !=
            GetFace(i+1,(Definitions::tri_index)(j+12)))
          std::cerr << "Outer R/L problem at : " << i <<
          ", " << j << std::endl;
        //test outer down
        if(i-w >= 0 &&
            GetFace(i,(Definitions::tri_index)(j+24)) !=
                GetFace(i-w,(Definitions::tri_index)(j+20)))
          std::cerr << "Outer D/U problem at : " << i <<
          ", " << j << std::endl;
        //test outer up
        if(GetFace(i,(Definitions::tri_index)(j+20)) !=
            GetFace(i+w,(Definitions::tri_index)(j+24)))
          std::cerr << "Outer U/D problem at : " << i <<
          ", " << j << std::endl;
        //test outer front
        if(i-w*h >= 0 &&
            GetFace(i,(Definitions::tri_index)(j+28)) !=
                GetFace(i-w*h,(Definitions::tri_index)(j+32)))
          std::cerr << "Outer F/B problem at : " << i <<
          ", " << j << std::endl;
        //test outer back
        if(GetFace(i,(Definitions::tri_index)(j+32)) !=
            GetFace(i+w*h,(Definitions::tri_index)(j+28)))
          std::cerr << "Outer B/F problem at : " << i <<
          ", " << j << std::endl;
        ////////////////////TESTING GET EDGES FROM FACE
        Definitions::tri_index f1,f2;
        std::array<CleaverCUDA::Edge*,3> e1, e2;
        //test outer left
        if(i >= 1) {
          f1 = (Definitions::tri_index)(j+12);
          f2 = (Definitions::tri_index)(j+16);
          e1 = GetFaceEdges(i,f1,(Definitions::tet_index)-1);
          e2 = GetFaceEdges(i-1,f2,(Definitions::tet_index)-1);
          for (size_t k = 0; k < 3; k++)
            if (e1[k] != e2[k])
              std::cerr << "Find edges problem left side at: " << i <<
              ", " << j << ", edge # " << k << std::endl;
        }
        //test outer right
        f1 = (Definitions::tri_index)(j+16);
        f2 = (Definitions::tri_index)(j+12);
        e1 = GetFaceEdges(i,f1,(Definitions::tet_index)-1);
        e2 = GetFaceEdges(i+1,f2,(Definitions::tet_index)-1);
        for (size_t k = 0; k < 3; k++)
          if (e1[k] != e2[k])
            std::cerr << "Find edges problem right side at: " << i <<
            ", " << j << ", edge # " << k << std::endl;
        //test outer down
        if(i >= w) {
          f1 = (Definitions::tri_index)(j+24);
          f2 = (Definitions::tri_index)(j+20);
          e1 = GetFaceEdges(i,f1,(Definitions::tet_index)-1);
          e2 = GetFaceEdges(i-w,f2,(Definitions::tet_index)-1);
          for (size_t k = 0; k < 3; k++)
            if (e1[k] != e2[k])
              std::cerr << "Find edges problem down side at: " << i <<
              ", " << j << ", edge # " << k << std::endl;
        }
        //test outer up
        f1 = (Definitions::tri_index)(j+20);
        f2 = (Definitions::tri_index)(j+24);
        e1 = GetFaceEdges(i,f1,(Definitions::tet_index)-1);
        e2 = GetFaceEdges(i+w,f2,(Definitions::tet_index)-1);
        for (size_t k = 0; k < 3; k++)
          if (e1[k] != e2[k])
            std::cerr << "Find edges problem up side at: " << i <<
            ", " << j << ", edge # " << k << std::endl;
        //test outer front
        if(i >= w*h) {
          f1 = (Definitions::tri_index)(j+28);
          f2 = (Definitions::tri_index)(j+32);
          e1 = GetFaceEdges(i,f1,(Definitions::tet_index)-1);
          e2 = GetFaceEdges(i-w*h,f2,(Definitions::tet_index)-1);
          for (size_t k = 0; k < 3; k++)
            if (e1[k] != e2[k])
              std::cerr << "Find edges problem front side at: " << i <<
              ", " << j << ", edge # " << k << std::endl;
        }
        //test outer back
        f1 = (Definitions::tri_index)(j+32);
        f2 = (Definitions::tri_index)(j+28);
        e1 = GetFaceEdges(i,f1,(Definitions::tet_index)-1);
        e2 = GetFaceEdges(i+w*h,f2,(Definitions::tet_index)-1);
        for (size_t k = 0; k < 3; k++)
          if (e1[k] != e2[k])
            std::cerr << "Find edges problem back side at: " << i <<
            ", " << j << ", edge # " << k << std::endl;
      }
    }
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
  } catch (const std::string& ex) {
    std::cerr << ex << std::endl;
  } catch (...) {
    std::cerr << "Unknown Exception" << std::endl;
  }
}
