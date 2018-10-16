/// \file ijkIO.txx
/// IO templates for reading/writing meshes.
/// - Input formats: Geomview .off.
/// - Output formats: Geomview .off, OpenInventor .iv (3D), Fig .fig (2D)
///   Stanford .ply, and Visualization Toolkit, .vtk.
/// - Version 0.1.5

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2018 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJKIO_
#define _IJKIO_

#include <iostream>
#include <vector>
#include <string>

#include "ijk.txx"

namespace IJK {

  // ******************************************
  // Miscellaneous local output routines
  // ******************************************

  // Forward declarations.
  namespace {

    template <typename CTYPE> void ijkoutVertexCoord
    (std::ostream & out, const int dim, const CTYPE * coord, const int numv);

    template <typename CTYPE, typename COLOR_TYPE> 
    void ijkoutVertexCoordColor
    (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

    template <typename VTYPE> void ijkoutPolygonVertices
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump);

    template <typename VTYPE> 
    void ijkoutPolygonVerticesRGB
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump,
     const unsigned char * poly_rgb);

    template <typename NTYPE0, typename NTYPE1,
              typename VTYPE, typename ITYPE>
    void ijkoutPolygonVerticesRGBII
    (std::ostream & out, 
     const NTYPE0 * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const NTYPE1 num_poly,
     const unsigned char * poly_rgb);

    template <typename NTYPE, typename VTYPE, typename ITYPE>
    void ijkoutPolytopeVertices
    (std::ostream & out, 
     const NTYPE * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const int num_poly);

    template <typename NTYPE, typename VTYPE, typename ITYPE,
              typename COLOR_TYPE>
    void ijkoutPolytopeVerticesColor
    (std::ostream & out, 
     const NTYPE * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const int num_poly,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

    template <typename VTYPE> void ijkoutQuadVertices
    (std::ostream & out, const VTYPE * quad_vert, const int numq,
     const bool flag_reorder_vertices);

    template <typename ITYPE, typename COLOR_TYPE>
    void ijkoutColor
    (std::ostream & out, const ITYPE i, 
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

    template <typename CTYPE0, typename CTYPE1, typename CTYPE2> 
    void ijkoutVectorEndpoints
    (std::ostream & out, const int dim, 
     const CTYPE0 * point_coord, const CTYPE1 * dir_coord, 
     const int num_vectors, const CTYPE2 scale);
  }

  // ******************************************
  // Write OpenInventor .iv file
  // ******************************************

  /// Output OpenInventor .iv file.
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * tri, const int numt)
    // out = output stream
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    IJK::PROCEDURE_ERROR error("ijkoutIV");

    if (dim != 3) {
      throw error
        ("Illegal dimension.  OpenInventor files are only for dimension 3.");
    }

    out << "#Inventor V2.1 ascii" << std::endl;
    out << std::endl;

    out << "Separator {" << std::endl;

    // set vertex ordering to clockwise to turn on two-sided lighting
    out << "  ShapeHints {" << std::endl;
    out << "    vertexOrdering CLOCKWISE" << std::endl;
    out << "  }" << std::endl;
    out << std::endl;

    out << "  IndexedFaceSet {" << std::endl;

    out << "    vertexProperty VertexProperty {" << std::endl;
    out << "      vertex [" << std::endl;

    // output vertex coordinates
    out << std::endl << "# vertex coordinates" << std::endl;
    for (int i = 0; i < numv; i++) {
      for (int d = 0; d < dim; d++) {
        out << coord[dim*i+d];
        if (d < dim-1) { out << " "; }
        else {	
          if (i < numv-1) { out << "," << std::endl; };
        };
      };
    };

    out << " ]" << std::endl;
    out << "    }" << std::endl;

    out << "    coordIndex [" << std::endl;
    // output triangle vertices
    out << std::endl << "# triangle vertices" << std::endl;
    for (int it = 0; it < numt; it++) {
      for (int d = 0; d < dim; d++) {
        out << tri[dim*it+d] << ",";
      };
      out << "-1";
      if (it < numt-1) {
        out << "," << std::endl;
      };
    };
    out << " ]" << std::endl;

    out << "  }" << std::endl;

    out << "}" << std::endl;
  }

  /// Output OpenInventor .iv file to standard output.
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (const int dim, const CTYPE * coord, const int numv,
   const VTYPE * tri, const int numt)
    // output OpenInventor .iv format to standard output
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    ijkoutIV(std::cout, dim, coord, numv, tri, numt); 
  }

  /// Output OpenInventor .iv file.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutIV(out, dim, vector2pointer(coord), coord.size()/dim,
             vector2pointer(simplex_vert), simplex_vert.size()/dim);
  }

  /// Output OpenInventor .iv file to standard output.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutIV(dim, vector2pointer(coord), coord.size()/dim,
             vector2pointer(simplex_vert), simplex_vert.size()/dim);
  }	

  // ******************************************
  // Write Geomview OFF file
  // ******************************************

  /// Output Geomview .off file.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param numv_per_simplex = Number of vertices per simplex.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices.
  ///        simplex_vert[numv_per_simplex*j+k] = k'th vertex of simplex j.
  /// @param nums = Number of simplices
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  {
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutPolygonVertices(out, numv_per_simplex, simplex_vert, nums);
  }

  /// Output Geomview .off file to standard output.
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, numv, 
              simplex_vert, nums);
  }

  /// \brief Output Geomview .off file.
  /// Every simplex has \a dim vertices.
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // set numv_per_simplex to dim
  {
    ijkoutOFF(out, dim, dim, coord, numv, simplex_vert, nums);
  }

  /// Output Geomview .off file to standard output.
  /// Every simplex has \a dim vertices.
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // output Geomview OFF format to standard output
    // set numv_per_simplex to dim
  {
    ijkoutOFF(dim, dim, coord, numv, simplex_vert, nums);
  }


  /// Output Geomview .off file.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF files
    // vector input format
  {
    ijkoutOFF
      (out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
       vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex);
  }

  /// Output Geomview .off file.
  /// Every simplex has \a dim vertices.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF files
    // vector input format
    // set numv_per_simplex to dim
  {
    int numv_per_simplex = dim;
    ijkoutOFF(out, dim, numv_per_simplex, coord, simplex_vert);
  }


  /// Output Geomview .off file to standard output.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  /// Output Geomview .off file to standard output.
  /// Every simplex has \a dim vertices.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF format to standard output
    // vector input format
    // set numv_per_simplex to dim
  {
    int numv_per_simplex = dim;
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  /// Output Geomview .off file.
  /// Two types of polytopes.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param poly1_vlist = List of vertices of poly 1.
  /// @param numv_per_poly1 = Number of vertices per polygon in poly1_vlist.
  /// @param num_poly1 = Number of poly1 polytopes.
  /// @param poly2_vlist = List of vertices of poly 2.
  /// @param numv_per_poly2 = Number of vertices per polygon in poly2_vlist.
  /// @param num_poly2 = Number of poly2 polytopes.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2> void ijkoutOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE1 * poly1_vlist, const int numv_per_poly1, const int num_poly1,
   const VTYPE2 * poly2_vlist, const int numv_per_poly2, const int num_poly2)
  {
    const int num_poly = num_poly1 + num_poly2;
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << num_poly << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutPolygonVertices(out, numv_per_poly1, poly1_vlist, num_poly1);
    ijkoutPolygonVertices(out, numv_per_poly2, poly2_vlist, num_poly2);
  }

  /// Output Geomview .off file.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2> void ijkoutOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE1> & poly1_vlist, const int numv_per_poly1,
   const std::vector<VTYPE2> & poly2_vlist, const int numv_per_poly2)
  {
    const int numv = coord.size()/dim;
    const int num_poly1 = poly1_vlist.size()/numv_per_poly1;
    const int num_poly2 = poly2_vlist.size()/numv_per_poly2;

    ijkoutOFF(out, dim, vector2pointer(coord), numv,
              vector2pointer(poly1_vlist), numv_per_poly1, num_poly1,
              vector2pointer(poly2_vlist), numv_per_poly2, num_poly2);
  }

  /// Output Geomview .off file.
  /// Different polytopes may have different numbers of vertices.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE> 
  void ijkoutPolytopeOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const NTYPE * num_poly_vert, const VTYPE * poly_vert,
   const ITYPE * first_poly_vert, const int num_poly)
  {
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << num_poly << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    ijkoutPolytopeVertices
      (out, num_poly_vert, poly_vert, first_poly_vert, num_poly);
  }

  /// Output Geomview .off file.
  /// - Different polytopes may have different numbers of vertices.
  /// - C++ STL vector format for coord[], num_poly_vert, poly_vert,
  ///   and first_poly_vert.
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE> 
  void ijkoutPolytopeOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<NTYPE> & num_poly_vert,
   const std::vector<VTYPE> & poly_vert,
   const std::vector<ITYPE> & first_poly_vert)
  {
    ijkoutPolytopeOFF
      (out, dim, vector2pointer(coord), coord.size()/dim,
       vector2pointer(num_poly_vert), vector2pointer(poly_vert),
       vector2pointer(first_poly_vert), num_poly_vert.size());
  }

  /// Output Geomview .off file. Color vertices header.
  inline void ijkoutColorVertOFFheader
  (std::ostream & out, const int dim,
   const int num_vert, const int num_poly,
   const bool is_back_colored)
  {
    out << "C";
    if (is_back_colored) out << "C";
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << num_vert << " " << num_poly << " " << 0 << std::endl;
  }

  /// Output Geomview .off file.
  /// Different polytopes may have different numbers of vertices.
  /// Color Vertices.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE,
            typename COLOR_TYPE> 
  void ijkoutPolyColorVertOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const NTYPE * num_poly_vert, const VTYPE * poly_vert,
   const ITYPE * first_poly_vert, const int num_poly,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPolyColorVertOFF");

    if (numv > 0 && front_color == NULL) {
      error.AddMessage("Programming error. front_color cannot be NULL.");
      throw error;
    }

    if (back_color == NULL)
      { ijkoutColorVertOFFheader(out, dim, numv, num_poly, false); }
    else
      { ijkoutColorVertOFFheader(out, dim, numv, num_poly, true); }


    ijkoutVertexCoordColor(out, dim, coord, numv, front_color, back_color);
    out << std::endl;

    ijkoutPolytopeVertices
      (out, num_poly_vert, poly_vert, first_poly_vert, num_poly);
  }

  /// Output Geomview .off file.
  /// - Different polytopes may have different numbers of vertices.
  /// - Color Vertices.
  /// - C++ STL vector format for num_poly_vert, poly_vert, first_poly_vert,
  ///        front_color[] and back_color[].
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE,
            typename COLOR_TYPE> 
  void ijkoutPolyColorVertOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const std::vector<NTYPE> & num_poly_vert, 
   const std::vector<VTYPE> & poly_vert,
   const std::vector<ITYPE> & first_poly_vert,
   const std::vector<COLOR_TYPE> & front_color, 
   const std::vector<COLOR_TYPE> & back_color)
  {
    ijkoutPolyColorVertOFF
      (out, dim, coord, numv, 
       vector2pointer(num_poly_vert), vector2pointer(poly_vert),
       vector2pointer(first_poly_vert), num_poly_vert.size(),
       vector2pointer(front_color), vector2pointer(back_color));
  }

  /// Output Geomview .off file. Color vertices.
  /// @param Output stream.
  /// @param dim Dimension.
  /// @param numv_per_simplex Number of vertices per simplex.
  /// @param coord[] coord[dim*i+k] = k'th coordinate of vertex j  (k < dim).
  /// @param numv Number of vertices.
  /// @param simplex_vert[] simplex_vert[dim*j+k] = 
  ///   k'th vertex index of simplex j.
  /// @param nums Number of simplices.
  /// @param front_color Array of front colors, 4 entries (RGBA) per vertex
  /// @parm back_color Array of backface colors, 4 entries (RGBA) per vertex
  ///                  May be NULL.
  template <typename T, typename COLOR_TYPE> void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    IJK::PROCEDURE_ERROR error("ijkoutColorVertOFF");

    if (numv > 0 && front_color == NULL) {
      error.AddMessage("Programming error. front_color cannot be NULL.");
      throw error;
    }

    if (back_color == NULL)
      { ijkoutColorVertOFFheader(out, dim, numv, nums, false); }
    else
      { ijkoutColorVertOFFheader(out, dim, numv, nums, true); }


    ijkoutVertexCoordColor(out, dim, coord, numv, front_color, back_color);
    ijkoutPolygonVertices(out, numv_per_simplex, simplex_vert, nums);
  }

  /// Output Geomview .off file to standard output. Color vertices.
  template <typename T, typename COLOR_TYPE> void ijkoutColorVertOFF
  (const int dim, const int numv_per_simplex,
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    ijkoutColorVertOFF
      (std::cout, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color vertices.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE, typename COLOR_TYPE> 
  void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    ijkoutColorVertOFF
      (out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
       vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color vertices.
  /// - C++ STL vector format for front_color[] and back_color[].
  template <typename T, typename COLOR_TYPE> void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const std::vector<COLOR_TYPE> & front_color, 
   const std::vector<COLOR_TYPE> & back_color)
  {
    ijkoutColorVertOFF
      (out, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       vector2pointer(front_color), vector2pointer(back_color));
  }

  /// Output Geomview .off file. Color vertices.
  /// Two types of polytopes.
  template <typename T, typename VTYPE1, typename VTYPE2,
            typename COLOR_TYPE> 
  void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const T * coord, const int numv,
   const VTYPE1 * poly1_vlist, const int numv_per_poly1, const int num_poly1,
   const VTYPE2 * poly2_vlist, const int numv_per_poly2, const int num_poly2,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    const int num_poly = num_poly1 + num_poly2;
    IJK::PROCEDURE_ERROR error("ijkoutColorVertOFF");

    if (numv > 0 && front_color == NULL) {
      error.AddMessage("Programming error. front_color cannot be NULL.");
      throw error;
    }

    if (back_color == NULL)
      { ijkoutColorVertOFFheader(out, dim, numv, num_poly, false); }
    else
      { ijkoutColorVertOFFheader(out, dim, numv, num_poly, true); }


    ijkoutVertexCoordColor(out, dim, coord, numv, front_color, back_color);
    ijkoutPolygonVertices(out, numv_per_poly1, poly1_vlist, num_poly1);
    ijkoutPolygonVertices(out, numv_per_poly2, poly2_vlist, num_poly2);
  }

  /// Output Geomview .off file.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2,
            typename COLOR_TYPE> 
  void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE1> & poly1_vlist, const int numv_per_poly1,
   const std::vector<VTYPE2> & poly2_vlist, const int numv_per_poly2,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    const int numv = coord.size()/dim;
    const int num_poly1 = poly1_vlist.size()/numv_per_poly1;
    const int num_poly2 = poly2_vlist.size()/numv_per_poly2;

    ijkoutColorVertOFF
      (out, dim, vector2pointer(coord), numv,
       vector2pointer(poly1_vlist), numv_per_poly1, num_poly1,
       vector2pointer(poly2_vlist), numv_per_poly2, num_poly2,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color faces header.
  inline void ijkoutColorFacesOFFheader
  (std::ostream & out, const int dim,
   const int num_vert, const int num_poly,
   const bool is_back_colored)
  {
    out << "D";
    if (is_back_colored) out << "D";
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << num_vert << " " << num_poly << " " << 0 << std::endl;
  }

  /// Output Geomview .off file.
  /// Different polytopes may have different numbers of vertices.
  /// Color faces.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE,
            typename COLOR_TYPE> 
  void ijkoutPolyColorFacesOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const NTYPE * num_poly_vert, const VTYPE * poly_vert,
   const ITYPE * first_poly_vert, const int num_poly,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPolyColorVertOFF");

    if (numv > 0 && front_color == NULL) {
      error.AddMessage("Programming error. front_color cannot be NULL.");
      throw error;
    }

    if (back_color == NULL)
      { ijkoutColorFacesOFFheader(out, dim, numv, num_poly, false); }
    else
      { ijkoutColorFacesOFFheader(out, dim, numv, num_poly, true); }


    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    ijkoutPolytopeVerticesColor
      (out, num_poly_vert, poly_vert, first_poly_vert, num_poly,
       front_color, back_color);
  }

  /// Output Geomview .off file.
  /// - Different polytopes may have different numbers of vertices.
  /// - Color faces.
  /// - C++ STL vector format for num_poly_vert, poly_vert, first_poly_vert,
  ///        front_color[] and back_color[].
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE,
            typename COLOR_TYPE> 
  void ijkoutPolyColorFacesOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const std::vector<NTYPE> & num_poly_vert, 
   const std::vector<VTYPE> & poly_vert,
   const std::vector<ITYPE> & first_poly_vert,
   const std::vector<COLOR_TYPE> & front_color, 
   const std::vector<COLOR_TYPE> & back_color)
  {
    ijkoutPolyColorFacesOFF
      (out, dim, coord, numv, 
       vector2pointer(num_poly_vert), vector2pointer(poly_vert),
       vector2pointer(first_poly_vert), num_poly_vert.size(),
       vector2pointer(front_color), vector2pointer(back_color));
  }

  /// Output Geomview .off file. Color simplices.
  /// @param Output stream.
  /// @param dim Dimension.
  /// @param numv_per_simplex Number of vertices per simplex.
  /// @param coord[] coord[dim*i+k] = k'th coordinate of vertex j  (k < dim).
  /// @param numv Number of vertices.
  /// @param simplex_vert[] simplex_vert[dim*j+k] = k'th vertex index of simplex j.
  /// @param nums Number of simplices.
  /// @param front_color Array of front colors, 4 entries (RGBA) per simplex.
  /// @parm back_color Array of backface colors, 4 entries (RGBA) per simplex.
  ///                  May be NULL.
  template <typename T, typename COLOR_TYPE> void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    IJK::PROCEDURE_ERROR error("ijkoutColorFacesOFF");

    if (nums > 0 && front_color == NULL) {
      error.AddMessage("Programming error. front_color cannot be NULL.");
      throw error;
    }

    if (back_color == NULL)
      { ijkoutColorFacesOFFheader(out, dim, numv, nums, false); }
    else
      { ijkoutColorFacesOFFheader(out, dim, numv, nums, true); }

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
        out << coord[iv*dim + d];
        if (d+1 < dim) out << " ";
      }
      out << std::endl;
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
        out << simplex_vert[is*numv_per_simplex + iv];
        if (iv < numv_per_simplex-1) { out << " "; }
      };
      out << "  ";
      for (int ic = 0; ic < 4; ic++) {
        out << front_color[4*is+ic];
        if (ic < 3) out << " ";
      }
      out << "  ";
      if (back_color != NULL) {
        for (int ic = 0; ic < 4; ic++) {
          out << back_color[4*is+ic];
          if (ic < 3) out << " ";
        }
      };
      out << std::endl;
    };

  }

  /// Output Geomview .off file. Color mesh faces.
  /// Two types of polytopes.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2,
            typename COLOR_TYPE> 
  void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE1 * poly1_vlist, const int numv_per_poly1, const int num_poly1,
   const COLOR_TYPE * poly1_front_color, const COLOR_TYPE * poly1_back_color,
   const VTYPE2 * poly2_vlist, const int numv_per_poly2, const int num_poly2,
   const COLOR_TYPE * poly2_front_color, const COLOR_TYPE * poly2_back_color)
  {
    const int num_poly = num_poly1 + num_poly2;
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    out << "D";
    if ((poly1_back_color != NULL) && (poly2_back_color != NULL))
      { out << "D"; }
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << num_poly << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutPolygonVerticesColor(out, numv_per_poly1, poly1_vlist, num_poly1, 
                               poly1_front_color, poly1_back_color);
    ijkoutPolygonVerticesColor(out, numv_per_poly2, poly2_vlist, num_poly2,
                               poly2_front_color, poly2_back_color);
  }

  /// Output Geomview .off file to standard output. Color simplices.
  template <typename T, typename COLOR_TYPE> void ijkoutColorFacesOFF
  (const int dim, const int numv_per_simplex, 
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    ijkoutColorFacesOFF
      (std::cout, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color simplices.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE, typename COLOR_TYPE> 
  void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
  {
    ijkoutColorFacesOFF
      (out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
       vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color simplices.
  /// - C++ STL vector format for front_color[] and back_color[].
  template <typename T, typename COLOR_TYPE> void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const std::vector<COLOR_TYPE> & front_color, 
   const std::vector<COLOR_TYPE> & back_color)
  {
    ijkoutColorFacesOFF
      (out, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       vector2pointer(front_color), vector2pointer(back_color));
  }

  /// Output Geomview .off file. Color faces.
  /// Two types of polytopes.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2, typename COLOR_TYPE> 
  void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE1> & poly1_vlist, const int numv_per_poly1,
   const COLOR_TYPE * poly1_front_color, const COLOR_TYPE * poly1_back_color,
   const std::vector<VTYPE2> & poly2_vlist, const int numv_per_poly2,
   const COLOR_TYPE * poly2_front_color, const COLOR_TYPE * poly2_back_color)
  {
    const int numv = coord.size()/dim;
    const int num_poly1 = poly1_vlist.size()/numv_per_poly1;
    const int num_poly2 = poly2_vlist.size()/numv_per_poly2;

    ijkoutColorFacesOFF
      (out, dim, vector2pointer(coord), numv,
       vector2pointer(poly1_vlist), numv_per_poly1, num_poly1,
       poly1_front_color, poly1_back_color,
       vector2pointer(poly2_vlist), numv_per_poly2, num_poly2,
       poly2_front_color, poly2_back_color);
  }

  /// Output Geomview .off file. Output vertex normals.
  template <typename CTYPE, typename NTYPE, typename VTYPE> 
  void ijkoutNormalsOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const CTYPE * coord, const NTYPE * normal, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // output Geomview OFF files with vertex normal information
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    if (dim == 3) { out << "NOFF" << std::endl; }
    else if (dim == 4) { out << "N4OFF" << std::endl;}
    else {
      out << "NnOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
        out << coord[iv*dim + d] << " ";
      }
      for (int d = 0; d < dim; d++) {
        out << normal[iv*dim+d];
        if (d < dim-1) { out << " "; }
        else { out << std::endl; };
      }
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
        out << simplex_vert[is*numv_per_simplex + iv];
        if (iv < numv_per_simplex-1) { out << " "; }
        else { out << std::endl; };
      };

    };

  }

  /// Output Geomview .off file. Output vertex normals.
  /// - C++ STL vector format for coord[], normal[] and simplex_vert[].
  template <typename CTYPE, typename NTYPE, typename VTYPE> 
  void ijkoutNormalsOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<NTYPE> & normal,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutNormalsOFF
      (out, dim, numv_per_simplex,
       vector2pointer(coord), vector2pointer(normal), coord.size()/dim,
       vector2pointer(simplex_vert),
       simplex_vert.size()/numv_per_simplex);
  }

  /// Output quadrilaterals to Geomview .off file.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param quad_vert = Array of quadrilateral vertices.
  ///        quad_vert[4*j+k] = k'th vertex of quad j.
  /// @param numq = Number of quadrilaterals.
  /// @param flag_reorder_vertices = Flag for reordering vertices.
  ///        If true, output vertices in counter-clockwise order around quad.
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (std::ostream & out, const int dim,
   const CTYPE * coord, const int numv,
   const VTYPE * quad_vert, const int numq,
   const bool flag_reorder_vertices)
  {
    IJK::PROCEDURE_ERROR error("ijkoutQuadOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << numq << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutQuadVertices
      (out, quad_vert, numq, flag_reorder_vertices);
  }

  /// Output quadrilaterals in Geomview .off format to standard output.
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (const int dim,
   const CTYPE * coord, const int numv,
   const VTYPE * quad_vert, const int numq,
   const bool flag_reorder_vertices)
  {
    ijkoutQuadOFF(std::cout, dim, coord, numv, quad_vert, numq);
  }

  /// Output quadrilaterals to Geomview .off
  /// - C++ STL vector format for coord[] and quad_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & quad_vert,
   const bool flag_reorder_vertices)
  {
    const int NUMV_PER_QUAD = 4;

    ijkoutQuadOFF(out, dim, vector2pointer(coord), coord.size()/dim,
                  vector2pointer(quad_vert), quad_vert.size()/NUMV_PER_QUAD,
                  flag_reorder_vertices);
  }

  /// Output quadrilaterals in Geomview .off format to standard output.
  /// - C++ STL vector format for coord[] and quad_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & quad_vert,
   const bool flag_reorder_vertices)
  {
    ijkoutQuadOFF(std::cout, dim, coord, quad_vert, flag_reorder_vertices);
  }

  // ******************************************
  // Write Geomview LINE file
  // ******************************************

  /// Output Geomview .line file.
  template <typename CTYPE, typename VTYPE> void ijkoutLINE
  (std::ostream & out, const int dim, 
   const CTYPE * coord, const int numv,
   const VTYPE * edge_vert, const int nume)
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // edge_vert[dim*j+k] = k'th vertex index of edge j
    // nums = number of simplices
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutLINE");

    if (dim == 3) { out << "LINE" << std::endl; }
    else if (dim == 4) { out << "4LINE" << std::endl;}
    else {
      error.AddMessage("Only dimensions 3 and 4 permitted for LINE output.");
      throw error;
    };

    out << numv << " " << nume << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    for (int ie = 0; ie < nume; ie++) {
      out << edge_vert[ie*NUM_EDGE_ENDPOINTS]
          << " " << edge_vert[ie*NUM_EDGE_ENDPOINTS+1] << std::endl;
    };

  }

  /// Output Geomview .line file.
  /// - C++ STL vector format for coord[] and edge_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutLINE
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & edge_vert)
  // output Geomview LINE file
  // vector input format
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutLINE");

    if (coord.size()%dim != 0) {
      error.AddMessage("Illegal number of vertex coordinates: ",
                       edge_vert.size(), ".");
      error.AddMessage
        ("  Number of coordinates must be divisible by the dimension: ",
         dim, ".");
      throw error;
    }

    if (edge_vert.size()%NUM_EDGE_ENDPOINTS != 0) {
      error.AddMessage("Illegal number of edge endpoints: ",
                       edge_vert.size(), ".");
      error.AddMessage("  Number of edge endpoints must be even.");
      throw error;
    }

    ijkoutLINE(out, dim, vector2pointer(coord), coord.size()/dim,
               vector2pointer(edge_vert), edge_vert.size()/NUM_EDGE_ENDPOINTS);
  }

  /// Output Geomview .line file color header
  template <typename COLOR_TYPE> 
  void ijkoutColorLINEheader
  (std::ostream & out, const int dim, const int numv, const int nume,
   const COLOR_TYPE rgba[4])
  {
    IJK::PROCEDURE_ERROR error("ijkoutColorLINE");

    if (dim == 3) { out << "LINEC" << std::endl; }
    else if (dim == 4) { out << "4LINEC" << std::endl;}
    else {
      error.AddMessage("Only dimensions 3 and 4 permitted for LINE output.");
      throw error;
    };

    out << rgba[0] << " " << rgba[1] << " " 
        << rgba[2] << " " << rgba[3] << std::endl;
    out << numv << " " << nume << std::endl;
  }


  /// Output Geomview .line file with color
  template <typename CTYPE, typename VTYPE, typename COLOR_TYPE> 
  void ijkoutColorLINE
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * edge_vert, const int nume, const COLOR_TYPE rgba[4])
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // edge_vert[dim*j+k] = k'th vertex index of edge j
    // nums = number of simplices
    // rgba[] = array of R,G,B,A values
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutColorLINE");

    ijkoutColorLINEheader(out, dim, numv, nume, rgba);

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    for (int ie = 0; ie < nume; ie++) {
      out << edge_vert[ie*NUM_EDGE_ENDPOINTS]
          << " " << edge_vert[ie*NUM_EDGE_ENDPOINTS+1] << std::endl;
    };

  }

  /// Output Geomview .line file.
  /// - C++ STL vector format for coord[] and edge_vert[].
  template <typename CTYPE, typename VTYPE, typename COLOR_TYPE> 
  void ijkoutColorLINE
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & edge_vert,
   const COLOR_TYPE rgba[4])
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutColorLINE");

    if (coord.size()%dim != 0) {
      error.AddMessage("Illegal number of vertex coordinates: ",
                       edge_vert.size(), ".");
      error.AddMessage
        ("  Number of coordinates must be divisible by the dimension: ",
         dim, ".");
      throw error;
    }

    if (edge_vert.size()%NUM_EDGE_ENDPOINTS != 0) {
      error.AddMessage("Illegal number of edge endpoints: ",
                       edge_vert.size(), ".");
      error.AddMessage("  Number of edge endpoints must be even.");
      throw error;
    }

    ijkoutColorLINE
      (out, dim, vector2pointer(coord), coord.size()/dim,
       vector2pointer(edge_vert), edge_vert.size()/NUM_EDGE_ENDPOINTS, rgba);
  }

  /// Output vectors to Geomview .line file.
  /// @param num_vectors Number of vectors
  /// @param scale Multiply all vectors by scale.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename COLOR_TYPE> 
  void ijkoutVectorsLINE
  (std::ostream & out, const int dim, 
   const CTYPE0 * point_coord, const CTYPE1 * dir_coord, 
   const int num_vectors, const CTYPE2 scale, const COLOR_TYPE rgba[4])
  {
    ijkoutColorLINEheader(out, dim, 2*num_vectors, num_vectors, rgba);

    ijkoutVectorEndpoints
      (out, dim, point_coord, dir_coord, num_vectors, scale);

    // Output line segments
    for (int iv = 0; iv < num_vectors; iv++) 
      { out << 2*iv << " " << 2*iv+1 << std::endl; }
  }

  /// Output vectors to Geomview .line file.
  /// @param num_vectors Number of vectors
  /// @param scale Multiply all vectors by scale.
  /// - C++ STL vector format for point_coord[] and dir_coord[].
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename COLOR_TYPE> 
  void ijkoutVectorsLINE
  (std::ostream & out, const int dim, 
   const std::vector<CTYPE0> & point_coord, 
   const std::vector<CTYPE1> & dir_coord, 
   const CTYPE2 scale, const COLOR_TYPE rgba[4])
  {
    IJK::PROCEDURE_ERROR error("ijkoutVectorsLINE");

    if (dim <= 0) {
      error.AddMessage("Illegal dimension: ", dim, ".");
      throw error;
    }

    if (point_coord.size() != dir_coord.size()) {
      error.AddMessage
        ("Error.  Number of points does not equal number of directions.");
      error.AddMessage("Number of points: ", point_coord.size()/dim, ".");
      error.AddMessage("Number of directions: ", dir_coord.size()/dim, ".");
      throw error;
    }

    ijkoutVectorsLINE
      (out, dim, vector2pointer(point_coord), vector2pointer(dir_coord),
       point_coord.size()/dim, scale, rgba);
  }

  // ******************************************
  // Read Geomview OFF file
  // ******************************************

  // local namespace
  namespace {
    inline void gobble_line(std::istream & in)
    {
      while (in.good() && !in.eof() && in.get() != '\n') {};
    }
  }

  /// \brief Read header of Geomview .off file.
  /// @param in Input stream.
  /// @param dim Vertex dimension.
  /// @param numv Number of vertices.
  /// @param nums Number of simplices.
  /// @param nume Number of edges.
  /// @param flag_normals Normals flag.  True if file contains vertex normals.
  inline void ijkinOFFheader
  (std::istream & in, int & dim, int & numv, int & nums, int & nume,
   bool & flag_normals)
  {
    std::string header_keyword;
    IJK::PROCEDURE_ERROR error("ijkinOFF");

    // Initialize
    dim = 0;
    numv = 0;
    nums = 0;
    nume = 0;
    flag_normals = false;

    if (!in.good()) {
      error.AddMessage("Error reading from input stream in.");
      throw error;
    }

    // read input
    in >> header_keyword;

    if (!in.good()) {
      error.AddMessage
        ("Error reading from header keyword from input stream in.");
      throw error;
    }

    bool flag_valid_header = true;
    int hksize = header_keyword.size() ;
    if (hksize < 3) { flag_valid_header = false; }
    else if (header_keyword.substr(hksize-3, 3) != "OFF")
      { flag_valid_header = false; }
    else {
      dim = 3;  // default dimension

      if (hksize > 3) {
        if (header_keyword.substr(hksize-4, 4) == "4OFF")
          { dim = 4; }
        else if (header_keyword.substr(hksize-4, 4) == "nOFF") 
          { in >> dim; }
      }

      std::string prefix = header_keyword.substr(0, hksize-3);

      int prefix_size = prefix.size();
      if ( prefix_size >= 1) {

        if (prefix[prefix_size-1] == '3' ||
            prefix[prefix_size-1] == '4' ||
            prefix[prefix_size-1] == 'n') {
          prefix = prefix.substr(0, prefix_size-1); 
          prefix_size = prefix_size-1;
        }
      }

      if (prefix_size >= 1) {

        if (prefix[prefix_size-1] == 'N') 
          { flag_normals = true; }
      }
    }

    if (!flag_valid_header) {
      error.AddMessage("Illegal Geomview .off file header: "+
                       header_keyword);
      throw error;
    };

    if (dim < 1)
      throw error("Dimension must be at least 1.");

    in >> numv;
    in >> nums;
    in >> nume;

    if (!in.good()) {
      error.AddMessage("Error reading number of vertices and polyhedra from input stream in.");
      throw error;
    }
  }

  /// \brief Read header of Geomview .off file.
  /// @param in = Input stream.
  /// @param dim = Vertex dimension.
  /// @param numv = Number of vertices.
  /// @param nums = Number of simplices.
  /// @param nume = Number of edges.
  inline void ijkinOFFheader
  (std::istream & in, int & dim, int & numv, int & nums, int & nume)
  {
    bool flag_normals;

    ijkinOFFheader(in, dim, numv, nums, nume, flag_normals);
  }

  /// \brief Read coordinates from Geomview .off file.
  ///
  /// Ignores any normal, color information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param numv = Number of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @pre Array coord[] is preallocated with size at least \a numv * \a dim.
  template <typename T> void ijkinOFFcoord
  (std::istream & in, const int dim, const int numv, T * coord)
  {
    IJK::PROCEDURE_ERROR error("ijkinOFFcoord");

    if (!in.good()) {
      error.AddMessage("Error reading coordinates from input stream in.");
      throw error;
    }

    if (coord == NULL) 
      if (numv > 0 && dim > 0) {
        error.AddMessage("Programming error. Array coord[] not allocated.");
        throw error;
      }

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) 
        { in >> coord[iv*dim + d]; };

      if (!in.good()) {
        error.AddMessage("Error reading coordinates of vertex ", iv,
                         " from input stream in.");
        throw error;
      }
    }
  }

  /// \brief Read coordinates and normals from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in Input stream.
  /// @param dim Dimension of vertices.
  /// @param numv Number of vertices.
  /// @param coord Array of coordinates. 
  ///        coord[dim*i+k] k'th coordinate of vertex i (k < dim).
  /// @pre Array coord[] is preallocated with size at least \a numv * \a dim.
  /// @param normal Array of normals.
  ///        normal[dim*i+k] k'th coordinate of normal of vertex i (k < dim).
  /// @pre Array normal[] is preallocated with size at least \a numv * \a dim.
  template <typename CTYPE, typename NTYPE> void ijkinNOFFcoord
  (std::istream & in, const int dim, const int numv, 
   CTYPE * coord, NTYPE * normal)
  {
    IJK::PROCEDURE_ERROR error("ijkinNOFFcoord");

    if (!in.good()) {
      error.AddMessage("Error reading coordinates from input stream in.");
      throw error;
    }

    if (coord == NULL) 
      if (numv > 0 && dim > 0) {
        error.AddMessage("Programming error. Array coord[] not allocated.");
        throw error;
      }

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) 
        { in >> coord[iv*dim + d]; };

      if (!in.good()) {
        error.AddMessage("Error reading coordinates of vertex ", iv,
                         " from input stream in.");
        throw error;
      }

      for (int d = 0; d < dim; d++) 
        { in >> normal[iv*dim + d]; };

      if (!in.good()) {
        error.AddMessage("Error reading normals of vertex ", iv,
                         " from input stream in.");
        throw error;
      }
    }
  }

  /// \brief Read coordinates from Geomview .off file.
  template <typename T> void ijkinOFFcoord
  (std::istream & in, const int dim, const int numv, std::vector<T> & coord)
  {
    if (dim < 1) {
      coord.clear();
      return;
    }

    coord.resize(dim*numv);

    ijkinOFFcoord(in, dim, numv, &(coord[0]));
  }

  /// \brief Read coordinates from Geomview .off file.
  template <typename CTYPE, typename NTYPE> void ijkinNOFFcoord
  (std::istream & in, const int dim, const int numv, 
   std::vector<CTYPE> & coord, std::vector<NTYPE> & normal)
  {
    if (dim < 1) {
      coord.clear();
      normal.clear();
      return;
    }

    coord.resize(dim*numv);
    normal.resize(dim*numv);

    ijkinNOFFcoord(in, dim, numv, &(coord[0]), &(normal[0]));
  }

  /// \brief Read vertex indices from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in = Input stream.
  /// @param numv = Number of vertices.
  /// @param vert = Array of vertices.
  ///        vert[k] = k'th vertex.
  /// @pre Array vert[] has been preallocated
  ///      with size at least \a numv.
  template <typename T> void ijkinOFFvert
  (std::istream & in, const int numv, T * vert)
  {
    for (int k = 0; k < numv; k++) 
      { in >> vert[k]; }
    gobble_line(in);
  }

  /// \brief Read polytope vertices of polytope \a js from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in = Input stream.
  /// @param jpoly = Polytope index.
  /// @param numv_per_poly = Number of vertices per polytope.
  /// @param poly_vert = Array of polytope vertices. 
  ///        poly_vert[numv_per_poly*js+k] = 
  ///            k'th vertex index of polytope jp.
  /// @pre Array poly_vert[] has been preallocated
  ///      with size at least \a numv_per_poly * \a (jpoly+1).
  template <typename T> void ijkinOFFpolyVert
  (std::istream & in, const int jpoly, const int numv_per_poly,
   T * poly_vert)
  {
    for (int k = 0; k < numv_per_poly; k++) 
      { in >> poly_vert[jpoly*numv_per_poly + k]; }
    gobble_line(in);
  }

  /// \brief Read simplex vertices of polytope \a jpoly from Geomview .off file.
  /// - C++ STL vector format for poly_vert[].
  /// @pre C++ vector poly_vert[] has size at least 
  ///        \a numv_per_poly * \a (jpoly+1).
  template <typename T> void ijkinOFFpolyVert
  (std::istream & in, const int jpoly, const int numv_per_poly,
   std::vector<T> & poly_vert)
  {
    ijkinOFFpolyVert(in, jpoly, numv_per_poly, &(poly_vert[0]));
  }

  /// \brief Read \a nums simplex vertices starting at simplex \a ifirst
  ///        from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in = Input stream.
  /// @param ifirst = Index of first simplex to read.
  /// @param nums = Number of simplices to read.
  /// @param numv_per_simplex = Number of vertices per simplex.
  /// @param simplex_vert = Array of simplex vertices. 
  ///        simplex_vert[numv_per_simplex*js+k] = 
  ///            k'th vertex index of simplex js.
  /// @pre Array simplex_vert[] has been preallocated
  ///      with size at least \a numv_per_simplex * \a (js+1).
  template <typename T> void ijkinOFFsimplex
  (std::istream & in, const int ifirst, const int nums, 
   const int numv_per_simplex, T * simplex_vert)
  {
    IJK::PROCEDURE_ERROR error("ijkinOFFsimplex");

    int nvert;
    for (int i = ifirst; i < ifirst + nums; i++) {

      in >> nvert;

      // Check that number of vertices is correct.
      if (nvert != numv_per_simplex) {
        error.AddMessage("Simplex ", i, " has ", nvert, " vertices.");
        error.AddMessage("All simplices are expected to have ",
                         numv_per_simplex, " vertices.");
        throw error;
      }

      ijkinOFFpolyVert(in, i, numv_per_simplex, simplex_vert);

      if (!in.good()) {
        error.AddMessage("Error reading vertices of simplex ", i, ".");
        throw error;
      }
    }
  }

  /// \brief Read \a nums simplex vertices starting at simplex \a ifirst
  ///        from Geomview .off file.
  /// - C++ STL vector format for simplex_vert[].
  /// @pre C++ vector simplex_vert[] has been preallocated
  ///      with size at least \a numv_per_simplex * \a (js+1).
  template <typename T> void ijkinOFFsimplex
  (std::istream & in, const int ifirst, const int nums, 
   const int numv_per_simplex, std::vector<T> & simplex_vert)
  {
    ijkinOFFsimplex(in, ifirst, nums, numv_per_simplex,
                    &(simplex_vert.front()));
  }


  /// \brief Read Geomview .off file.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param mesh_dim = Dimension of mesh.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices. 
  ///        simplex_vert[dim*j+k] = k'th vertex index of simplex j.
  /// @param nums = Number of simplices.
  /// @pre All simplices have the same dimension.
  template <typename T> void ijkinOFF
  (std::istream & in, int & dim, int & mesh_dim,
   T * & coord, int & numv, int * & simplex_vert, int & nums)
  {
    int nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord = NULL;
    simplex_vert = NULL;
    mesh_dim = 0;            // default mesh dimension

    ijkinOFFheader(in, dim, numv, nums, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    if (nums > 0) {
      // use first simplex to set mesh dimension
      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert > 0) { 
        mesh_dim = num_simplex_vert-1; 
      }
      else { mesh_dim = 0; };

      simplex_vert = new int[nums*num_simplex_vert];

      // read in first simplex
      ijkinOFFpolyVert(in, 0, num_simplex_vert, simplex_vert);

      ijkinOFFsimplex(in, 1, nums-1, num_simplex_vert, simplex_vert);
    };
  }

  /// \brief Read Geomview .off file.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename T> void ijkinOFF
  (std::istream & in, int & dim, int & mesh_dim,
   std::vector<T> & coord, std::vector<int> & simplex_vert)
  {
    int numv, nums, nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord.clear();
    simplex_vert.clear();

    mesh_dim = 0;            // default mesh dimension

    ijkinOFFheader(in, dim, numv, nums, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    if (nums > 0) {
      // use first simplex to set mesh dimension
      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert > 0) { 
        mesh_dim = num_simplex_vert-1; 
      }
      else { mesh_dim = 0; };

      simplex_vert.resize(nums*num_simplex_vert);

      // read vertices of first simplex
      ijkinOFFpolyVert(in, 0, num_simplex_vert, simplex_vert);

      // read remaining simplex vertices
      ijkinOFFsimplex(in, 1, nums-1, num_simplex_vert, simplex_vert);
    };
  }

  /// Read Geomview .off file from standard input.
  template <typename T> void ijkinOFF
  (int & dim, int & mesh_dim,
   T * & coord, int & numv, int * & simplex_vert, int & nums)
    // input Geomview OFF format to standard input
    // in = input stream
    // dim = dimension
    // mesh_dim = mesh dimension
    //   Precondition: All simplices have the same dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    ijkinOFF(std::cin, dim, mesh_dim, coord, numv, simplex_vert, nums);
  }


  /// \brief Read Geomview .off file where each simplex has dimension \a dim-1.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices. 
  ///        simplex_vert[dim*j+k] = k'th vertex index of simplex j.
  /// @param nums = Number of simplices.
  /// @pre All simplices have the same dimension.
  template <typename T> void ijkinOFF
  (std::istream & in, int & dim, T * & coord, int & numv,
   int * & simplex_vert, int & nums)
  {
    int nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord = NULL;
    simplex_vert = NULL;

    if (!in.good()) {
      throw error("Error: Corrupted input stream. Unable to read input.");
    }

    ijkinOFFheader(in, dim, numv, nums, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    simplex_vert = new int[nums*dim];
    int numv_per_simplex = dim;
    ijkinOFFsimplex(in, 0, nums, numv_per_simplex, simplex_vert);
  }

  /// \brief Read Geomview .off file where each simplex has dimension \a dim-1 
  ///   from standard input.
  template <typename T> void ijkinOFF
  (int & dim, T * & coord, int & numv, int * & simplex_vert, int & nums)
  {
    ijkinOFF(std::cin, dim, coord, numv, simplex_vert, nums);
  }

  /// \brief Read Geomview .off file.
  /// - Read normal information from files with NOFF header.
  /// - C++ STL vector format for coord[], normal[], and simplex_vert[].
  template <typename CTYPE, typename NTYPE> void ijkinOFF
  (std::istream & in, int & dim, int & mesh_dim,
   std::vector<CTYPE> & coord, std::vector<NTYPE> & normal, 
   std::vector<int> & simplex_vert)
  {
    bool flag_normals;
    int numv, nums, nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord.clear();
    normal.clear();
    simplex_vert.clear();

    mesh_dim = 0;            // default mesh dimension

    ijkinOFFheader(in, dim, numv, nums, nume, flag_normals);

    if (flag_normals) {
      ijkinNOFFcoord(in, dim, numv, coord, normal);
    }
    else {
      ijkinOFFcoord(in, dim, numv, coord);
    }

    if (nums > 0) {
      // use first simplex to set mesh dimension
      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert > 0) { 
        mesh_dim = num_simplex_vert-1; 
      }
      else { mesh_dim = 0; };

      simplex_vert.resize(nums*num_simplex_vert);

      // read vertices of first simplex
      ijkinOFFpolyVert(in, 0, num_simplex_vert, simplex_vert);

      // read remaining simplex vertices
      ijkinOFFsimplex(in, 1, nums-1, num_simplex_vert, simplex_vert);
    };
  }

  /// \brief Read triangles and quadrilaterals from Geomview .off file.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param tri_vert = Array of triangle vertices. 
  ///        tri_vert[3*j+k] = k'th vertex index of triangle j.
  /// @param nums = Number of simplices.
  /// @param quad_vert = Array of quadrilateral vertices.
  ///        simplex_vert[4*j+k] = k'th vertex index of quad j.
  /// @param numq = Number of quadrilaterals.
  template <typename T> 
  void ijkinOFF_tri_quad
  (std::istream & in, int & dim, std::vector<T> & coord, 
   std::vector<int> & tri_vert, std::vector<int> & quad_vert)
  {
    const int NUM_TRIANGLE_VERTICES = 3;
    const int NUM_QUAD_VERTICES = 4;
    IJK::PROCEDURE_ERROR error("ijkinOFF_tri_quad");

    int nume;

    coord.clear();
    tri_vert.clear();
    quad_vert.clear();

    // nump: Total number of triangle and quads.
    int numv, nump;
    ijkinOFFheader(in, dim, numv, nump, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    int numt = 0;
    int numq = 0;
    for (int i = 0; i < nump; i++) {
      int nvert;

      in >> nvert;

      if (nvert == NUM_TRIANGLE_VERTICES) {
        tri_vert.resize(NUM_TRIANGLE_VERTICES*(numt+1));
        ijkinOFFpolyVert(in, numt, NUM_TRIANGLE_VERTICES, tri_vert);
        numt++;
      }
      else if (nvert == NUM_QUAD_VERTICES) {
        quad_vert.resize(NUM_QUAD_VERTICES*(numq+1));
        ijkinOFFpolyVert(in, numq, NUM_QUAD_VERTICES, quad_vert);
        numq++;
      }
      else {
        // Illegal number of vertices.
        error.AddMessage("Polygon ", i, " has ", nvert, " vertices.");
        error.AddMessage("All polygons are expected to have 3 or 4 vertices.");
        throw error;
      }

    }

  }

  /// \brief Read polygons from Geomview .off file.
  /// \brief Read triangles and quadrilaterals into separate arrays.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param tri_vert = Array of triangle vertices. 
  ///        tri_vert[3*j+k] = k'th vertex index of triangle j.
  /// @param quad_vert = Array of quadrilateral vertices.
  ///        quad_vert[4*j+k] = k'th vertex index of quad j.
  /// @param num_poly_vert = Array of number of polygon vertices.
  ///        num_poly_vert[i] = Number of vertices of polygon i.
  /// @param poly_vert = Array of polygon vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polygon.
  ///        Polygon i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename T> 
  void ijkinPolyOFF
  (std::istream & in, int & dim, std::vector<T> & coord, 
   std::vector<int> & tri_vert, std::vector<int> & quad_vert,
   std::vector<int> & num_poly_vert,
   std::vector<int> & poly_vert, std::vector<int> & first_poly_vert)
  {
    const int NUM_TRIANGLE_VERTICES = 3;
    const int NUM_QUAD_VERTICES = 4;
    IJK::PROCEDURE_ERROR error("ijkinPolyOFF");

    int nume;

    coord.clear();
    tri_vert.clear();
    quad_vert.clear();
    num_poly_vert.clear();
    poly_vert.clear();
    first_poly_vert.clear();

    // nump: Total number of polygons.
    int numv, nump;
    ijkinOFFheader(in, dim, numv, nump, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    int numt = 0;
    int numq = 0;
    int num_not_qt = 0;   // Number of polygons which are not tri or quad
    for (int i = 0; i < nump; i++) {
      int nvert;

      in >> nvert;

      if (nvert == NUM_TRIANGLE_VERTICES) {
        tri_vert.resize(NUM_TRIANGLE_VERTICES*(numt+1));
        ijkinOFFpolyVert(in, numt, NUM_TRIANGLE_VERTICES, tri_vert);
        numt++;
      }
      else if (nvert == NUM_QUAD_VERTICES) {
        quad_vert.resize(NUM_QUAD_VERTICES*(numq+1));
        ijkinOFFpolyVert(in, numq, NUM_QUAD_VERTICES, quad_vert);
        numq++;
      }
      else {
        num_poly_vert.push_back(nvert);
        int k = poly_vert.size();
        first_poly_vert.push_back(k);
        poly_vert.resize(k+nvert);
        ijkinOFFvert(in, nvert, &(poly_vert[k]));
        num_not_qt++;
      }
    }
  }

  /// \brief Read list of polytopes from Geomview .off file.
  /// @param in = Input stream.
  /// @param num_poly_vert = Array of number of polytope vertices.
  ///        num_poly_vert[i] = Number of vertices of polytope i.
  /// @param poly_vert = Array of polytope vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polytope.
  ///        Polytope i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename NTYPE, typename VTYPE, typename ITYPE>
  void ijkinOFFpolyList
  (std::istream & in, const int nump, std::vector<NTYPE> & num_poly_vert,
   std::vector<VTYPE> & poly_vert, std::vector<ITYPE> & first_poly_vert)
  {
    for (int i = 0; i < nump; i++) {
      NTYPE nvert;

      in >> nvert;
      num_poly_vert.push_back(nvert);
      int k = poly_vert.size();
      first_poly_vert.push_back(k);
      poly_vert.resize(k+nvert);
      ijkinOFFvert(in, nvert, &(poly_vert[k]));
    }
  }

  /// \brief Read polytopes from Geomview .off file.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param num_poly_vert = Array of number of polytope vertices.
  ///        num_poly_vert[i] = Number of vertices of polytope i.
  /// @param poly_vert = Array of polytope vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polytope.
  ///        Polytope i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename T> 
  void ijkinPolytopeOFF
  (std::istream & in, int & dim, T * & coord, int & numv, 
   std::vector<int> & num_poly_vert,
   std::vector<int> & poly_vert, std::vector<int> & first_poly_vert)
  {
    IJK::PROCEDURE_ERROR error("ijkinPolytopeOFF");

    int nume, nump;

    coord = NULL;
    num_poly_vert.clear();
    poly_vert.clear();
    first_poly_vert.clear();

    ijkinOFFheader(in, dim, numv, nump, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    ijkinOFFpolyList
      (in, nump, num_poly_vert, poly_vert, first_poly_vert);
  }

  /// \brief Read polytopes from Geomview .off file.
  /// - C++ STL vector format for coord[].
  /// - Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param num_poly_vert = Array of number of polytope vertices.
  ///        num_poly_vert[i] = Number of vertices of polytope i.
  /// @param poly_vert = Array of polytope vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polytope.
  ///        Polytope i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename T> 
  void ijkinPolytopeOFF
  (std::istream & in, int & dim, std::vector<T> & coord, 
   std::vector<int> & num_poly_vert,
   std::vector<int> & poly_vert, std::vector<int> & first_poly_vert)
  {
    IJK::PROCEDURE_ERROR error("ijkinPolytopeOFF");

    int nume;

    coord.clear();
    num_poly_vert.clear();
    poly_vert.clear();
    first_poly_vert.clear();

    // nump: Total number of polytopes.
    int numv, nump;
    ijkinOFFheader(in, dim, numv, nump, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    ijkinOFFpolyList
      (in, nump, num_poly_vert, poly_vert, first_poly_vert);
  }


  // ******************************************
  // Read Geomview LINE file
  // ******************************************

  /// \brief Read header of Geomview .line file.
  /// Ignores line color.
  /// @param in = Input stream.
  /// @param dim = Vertex dimension.
  /// @param numv = Number of vertices.
  /// @param nume = Number of edges.
  inline void ijkinLINEheader
  (std::istream & in, int & dim, int & numv, int & nume)
  {
    float r,g,b,a;
    std::string header_keyword;
    IJK::PROCEDURE_ERROR error("ijkinLINE");
    if (!in.good()) {
      error.AddMessage("Error reading from input stream in.");
      throw error;
    }

    // read input
    in >> header_keyword;

    if (!in.good()) {
      error.AddMessage
        ("Error reading from header keyword from input stream in.");
      throw error;
    }

    if (header_keyword == "LINEC") {
      dim = 3;
      in >> r >> g >> b >> a;

      if (!in.good()) {
        error.AddMessage("Error reading number of rgba edge color from input stream in.");
        throw error;
      }
    }
    else if (header_keyword == "LINE") {
      dim = 3;
    }
    else {
      error.AddMessage("Illegal Geomview .off file header: "+
                       header_keyword);
      throw error;
    };

    if (dim < 1)
      throw error("Dimension must be at least 1.");

    in >> numv;
    in >> nume;

    if (!in.good()) {
      error.AddMessage("Error reading number of vertices and edges from input stream in.");
      throw error;
    }
  }

  /// \brief Read \a nume edges from Geomview .off file.
  ///
  /// @param in = Input stream.
  /// @param nume = Number of edges to read.
  /// @param edge_endpoint = Array of edge endpoints. 
  ///        edge_enpdoint[2*je+k] = index of k'th endpoint of edge je.
  /// @pre Array edge_endpoint[] has been preallocated
  ///      with size at least 2 * \a nume.
  template <typename T> void ijkinLINEedge
  (std::istream & in, const int nume, T * & edge_endpoint)
  {
    IJK::PROCEDURE_ERROR error("ijkinLINEedge");

    for (int i = 0; i < nume; i++) {

      for (int k = 0; k < 2; k++) 
        { in >> edge_endpoint[2*i + k]; }
      gobble_line(in);

      if (!in.good()) {
        error.AddMessage("Error reading endpoints of edge ", i, ".");
        throw error;
      }
    }
  }

  /// \brief Read Geomview .line file.
  ///
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param edge_endpoint = Array of edge endpoints. 
  ///        edge_endpoint[dim*j+k] = Vertex index of endpoint k of edge j.
  /// @param nume = Number of edges.
  template <typename T> void ijkinLINE
  (std::istream & in, int & dim, T * & coord, int & numv, 
   int * & edge_endpoint, int & nume)
  {
    IJK::PROCEDURE_ERROR error("ijkinLINE");

    coord = NULL;
    edge_endpoint = NULL;

    ijkinLINEheader(in, dim, numv, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    int num_edge_endpoint = 0;
    edge_endpoint = new int[2*nume];

    // read in first edge
    ijkinLINEedge(in, nume, edge_endpoint);
  }


  // ******************************************
  // Write .ply file
  // ******************************************

  /// Output .ply file header
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param numv = Number of vertices.
  /// @param numf = Number of faces.
  inline void ijkoutPLYheader
  (std::ostream & out, const int dim, const int numv, const int numf)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPLYheader");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex " << numv << std::endl;
    out << "property float x" << std::endl;
    out << "property float y" << std::endl;
    out << "property float z" << std::endl;
    out << "element face " << numf << std::endl;
    out << "property list uchar int vertex_index" << std::endl;
    out << "end_header" << std::endl;
  }


  /// Output .ply file.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param numv_per_simplex = Number of vertices per simplex.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices.
  ///        simplex_vert[numv_per_simplex*j+k] = k'th vertex of simplex j.
  /// @param nums = Number of simplices
  template <typename CTYPE, typename VTYPE> void ijkoutPLY
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPLY");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    ijkoutPLYheader(out, dim, numv, nums);
    ijkoutVertexCoord(out, dim, coord, numv);
    ijkoutPolygonVertices(out, numv_per_simplex, simplex_vert, nums);
  }


  /// Output .ply file.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutPLY
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutPLY
      (out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
       vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex);
  }


  /// Output .ply file to standard output.
  template <typename CTYPE, typename VTYPE> void ijkoutPLY
  (const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  {
    ijkoutPLY(std::cout, dim, numv_per_simplex, coord, numv,
              simplex_vert, nums);
  }


  /// Output .ply file to standard output.
  /// - C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutPLY
  (const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutPLY(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  /// Output .ply file.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param numv_per_simplex = Number of vertices per simplex.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param quad_vert = Array of quadrilateral vertices.
  ///        quad_vert[4*j+k] = k'th vertex of quad j.
  /// @param numq = Number of quadrilaterals.
  /// @param flag_reorder_vertices = Flag for reordering vertices.
  ///        If true, output vertices in counter-clockwise order around quad.
  template <typename CTYPE, typename VTYPE> void ijkoutQuadPLY
  (std::ostream & out, const int dim,
   const CTYPE * coord, const int numv,
   const VTYPE * quad_vert, const int numq,
   const bool flag_reorder_vertices)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPLY");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    ijkoutPLYheader(out, dim, numv, numq);
    ijkoutVertexCoord(out, dim, coord, numv);
    ijkoutQuadVertices
      (out, quad_vert, numq, flag_reorder_vertices);
  }

  /// Output quadrilaterals to .ply file
  /// - C++ STL vector format for coord[] and quad_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutQuadPLY
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & quad_vert,
   const bool flag_reorder_vertices)
  {
    const int NUMV_PER_QUAD = 4;

    ijkoutQuadPLY(out, dim, vector2pointer(coord), coord.size()/dim,
                  vector2pointer(quad_vert), quad_vert.size()/NUMV_PER_QUAD,
                  flag_reorder_vertices);
  }

  /// Output .ply file.
  /// Two types of polytopes.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param poly1_vlist = List of vertices of poly 1.
  /// @param numv_per_poly1 = Number of vertices per polygon in poly1_vlist.
  /// @param num_poly1 = Number of poly1 polytopes.
  /// @param poly2_vlist = List of vertices of poly 2.
  /// @param numv_per_poly2 = Number of vertices per polygon in poly2_vlist.
  /// @param num_poly2 = Number of poly2 polytopes.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2> void ijkoutPLY
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE1 * poly1_vlist, const int numv_per_poly1, const int num_poly1,
   const VTYPE2 * poly2_vlist, const int numv_per_poly2, const int num_poly2)
  {
    const int num_poly = num_poly1 + num_poly2;
    IJK::PROCEDURE_ERROR error("ijkoutPLY");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    ijkoutPLYheader(out, dim, numv, num_poly);
    ijkoutVertexCoord(out, dim, coord, numv);
    ijkoutPolygonVertices(out, numv_per_poly1, poly1_vlist, num_poly1);
    ijkoutPolygonVertices(out, numv_per_poly2, poly2_vlist, num_poly2);
  }

  /// Output .ply file.
  /// - C++ STL vector format for coord[], poly1_vlist[] and poly2_vlist[].
  template <typename CTYPE, typename VTYPE1, typename VTYPE2> void ijkoutPLY
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE1> & poly1_vlist, const int numv_per_poly1,
   const std::vector<VTYPE2> & poly2_vlist, const int numv_per_poly2)
  {
    const int numv = coord.size()/dim;
    const int num_poly1 = poly1_vlist.size()/numv_per_poly1;
    const int num_poly2 = poly2_vlist.size()/numv_per_poly2;

    ijkoutPLY(out, dim, vector2pointer(coord), numv,
              vector2pointer(poly1_vlist), numv_per_poly1, num_poly1,
              vector2pointer(poly2_vlist), numv_per_poly2, num_poly2);
  }


  /// Output .ply file.
  /// Different polytopes may have different numbers of vertices.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE> 
  void ijkoutPolytopePLY
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const NTYPE * num_poly_vert, const VTYPE * poly_vert,
   const ITYPE * first_poly_vert, const int num_poly)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPolytopePLY");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    ijkoutPLYheader(out, dim, numv, num_poly);
    ijkoutVertexCoord(out, dim, coord, numv);
    ijkoutPolytopeVertices
      (out, num_poly_vert, poly_vert, first_poly_vert, num_poly);
  }


  // ******************************************
  // Write .ply file with colored objects
  // ******************************************

  /// Output .ply file header with faces colored.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param numv = Number of vertices.
  /// @param nums = Number of simplices
  inline void ijkoutColorFacesPLYheader
  (std::ostream & out, const int dim, const int numv, const int nums)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPLYheader");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex " << numv << std::endl;
    out << "property float x" << std::endl;
    out << "property float y" << std::endl;
    out << "property float z" << std::endl;
    out << "element face " << nums << std::endl;
    out << "property list uchar int vertex_index" << std::endl;
    out << "property uchar red" << std::endl;
    out << "property uchar green" << std::endl;
    out << "property uchar blue" << std::endl;
    out << "end_header" << std::endl;
  }


  /// Output .ply file with colored faces.
  /// Two types of polytopes.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param poly1_vlist = List of vertices of poly 1.
  /// @param numv_per_poly1 = Number of vertices per polygon in poly1_vlist.
  /// @param num_poly1 = Number of poly1 polytopes.
  /// @param poly1_rgb[] Must be from 0 to 255.
  /// @param poly2_vlist = List of vertices of poly 2.
  /// @param numv_per_poly2 = Number of vertices per polygon in poly2_vlist.
  /// @param num_poly2 = Number of poly2 polytopes.
  /// @param poly2_rgb[] Must be from 0 to 255.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2>
  void ijkoutColorFacesPLY
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE1 * poly1_vlist, const int numv_per_poly1, const int num_poly1,
   const unsigned char * poly1_rgb,
   const VTYPE2 * poly2_vlist, const int numv_per_poly2, const int num_poly2,
   const unsigned char * poly2_rgb)
  {
    const int num_poly = num_poly1 + num_poly2;
    IJK::PROCEDURE_ERROR error("ijkoutColorFacesPLY");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    if (num_poly1 > 0 && poly1_vlist == NULL) {
      error.AddMessage
        ("Programming error.  Number of poly1 is ", num_poly1,
         " but poly1 vertex list is empty.");
      throw error;
    }

    if (num_poly2 > 0 && poly2_vlist == NULL) {
      error.AddMessage
        ("Programming error.  Number of poly2 is ", num_poly2,
         " but poly2 vertex list is empty.");
      throw error;
    }

    ijkoutColorFacesPLYheader(out, dim, numv, num_poly);
    ijkoutVertexCoord(out, dim, coord, numv);
    ijkoutPolygonVerticesRGB(out, numv_per_poly1, poly1_vlist, num_poly1, 
                             poly1_rgb);
    ijkoutPolygonVerticesRGB(out, numv_per_poly2, poly2_vlist, num_poly2, 
                             poly2_rgb);
  }


  /// Output .ply file with colored faces.
  /// - Two types of polytopes.
  /// - C++ STL vector format for coord[], poly1_vlist[], poly2_vlist[]
  ///     poly1_color[] and poly2_color[].
  template <typename CTYPE, typename VTYPE1, typename VTYPE2>
  void ijkoutColorFacesPLY
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE1> & poly1_vlist, const int numv_per_poly1,
   const std::vector<unsigned char> poly1_color,
   const std::vector<VTYPE2> & poly2_vlist, const int numv_per_poly2,
   const std::vector<unsigned char> poly2_color)
  {
    const int numv = coord.size()/dim;
    const int num_poly1 = poly1_vlist.size()/numv_per_poly1;
    const int num_poly2 = poly2_vlist.size()/numv_per_poly2;

    ijkoutColorFacesPLY
      (out, dim, vector2pointer(coord), numv,
       vector2pointer(poly1_vlist), numv_per_poly1, num_poly1,
       vector2pointer(poly1_color),
       vector2pointer(poly2_vlist), numv_per_poly2, num_poly2,
       vector2pointer(poly2_color));
  }


  /// Output .ply file with colored faces.
  /// - Version where polygons may have different numbers of vertices.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE> 
  void ijkoutColorFacesPLYII
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const NTYPE * num_poly_vert, const VTYPE * poly_vert,
   const ITYPE * first_poly_vert, const int num_poly,
   const unsigned char * poly_rgb)
  {
    IJK::PROCEDURE_ERROR error("ijkoutColorFacesPLYII");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    ijkoutColorFacesPLYheader(out, dim, numv, num_poly);
    ijkoutVertexCoord(out, dim, coord, numv);
    ijkoutPolygonVerticesRGBII
      (out, num_poly_vert, poly_vert, first_poly_vert, num_poly, poly_rgb);

  }


  // ******************************************
  // Write .ply edge file
  // ******************************************

  /// Output .ply file header of edges
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param numv = Number of vertices.
  /// @param nums = Number of simplices
  inline void ijkoutEdgePLYheader
  (std::ostream & out, const int dim, const int numv, const int nume)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPLYheader");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex " << numv << std::endl;
    out << "property float x" << std::endl;
    out << "property float y" << std::endl;
    out << "property float z" << std::endl;
    out << "element edge " << nume << std::endl;
    out << "property int vertex1" << std::endl;
    out << "property int vertex2" << std::endl;
    out << "end_header" << std::endl;
  }

  /// Output .ply files of edges.
  /// out = output stream
  /// dim = dimension
  /// coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
  /// numv = number of vertices 
  /// edge_vert[dim*j+k] = k'th vertex index of edge j
  /// nume = number of edges
  template <typename CTYPE, typename VTYPE> void ijkoutEdgePLY
  (std::ostream & out, const int dim, 
   const CTYPE * coord, const int numv,
   const VTYPE * edge_vert, const int nume)
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutEdgePLY");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 implemented for .ply files.");
      throw error;
    }

    ijkoutEdgePLYheader(out, dim, numv, nume);
    ijkoutVertexCoord(out, dim, coord, numv);

    for (int ie = 0; ie < nume; ie++) {
      out << edge_vert[ie*NUM_EDGE_ENDPOINTS]
          << " " << edge_vert[ie*NUM_EDGE_ENDPOINTS+1] << std::endl;
    };

  }


  /// Output .ply edge file.
  /// - C++ STL vector format for coord[] and edge_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutEdgePLY
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & edge_vert)
  // output Geomview LINE file
  // vector input format
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutEdgePLY");

    if (coord.size()%dim != 0) {
      error.AddMessage("Illegal number of vertex coordinates: ",
                       edge_vert.size(), ".");
      error.AddMessage
        ("  Number of coordinates must be divisible by the dimension: ",
         dim, ".");
      throw error;
    }

    if (edge_vert.size()%NUM_EDGE_ENDPOINTS != 0) {
      error.AddMessage("Illegal number of edge endpoints: ",
                       edge_vert.size(), ".");
      error.AddMessage("  Number of edge endpoints must be even.");
      throw error;
    }

    ijkoutEdgePLY
      (out, dim, vector2pointer(coord), coord.size()/dim,
       vector2pointer(edge_vert), edge_vert.size()/NUM_EDGE_ENDPOINTS);
  }


  // ******************************************
  // Write .vtk file
  // ******************************************

  /// Output ASCII .vtk file header for unstructured grid.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  inline void ijkoutVTKheader
  (std::ostream & out, const char * dataset_name, const int dim)
  {
    IJK::PROCEDURE_ERROR error("ijkoutVTKheader");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 available for .vtk files.");
      throw error;
    }

    using std::endl;

    out << "# vtk DataFile Version 4.1" << endl;
    out << dataset_name << endl;
    out << "ASCII" << endl;
  }


  /// Output cell types.
  inline void ijkoutCellTypesVTK
  (std::ostream & out, const int num_cells, const int itype)
  {
    using std::endl;

    out << "CELL_TYPES " << num_cells << endl;
    for (int i = 0; i < num_cells; i++) {
      out << itype << endl;
    }
  }


  /// Output hexahedra types.
  inline void ijkoutHexahedraTypesVTK
  (std::ostream & out, const int num_hex)
  {
    const int HEXAHEDRON_TYPE(12);

    ijkoutCellTypesVTK(out, num_hex, HEXAHEDRON_TYPE);
  }


  /// Output tetrahedra types.
  inline void ijkoutTetrahedraTypesVTK
  (std::ostream & out, const int num_hex)
  {
    const int TETRAHEDRON_TYPE(10);

    ijkoutCellTypesVTK(out, num_hex, TETRAHEDRON_TYPE);
  }


  /// Output hex vertices, reordering in VTK order.
  template <typename VTYPE, typename NTYPE> void ijkoutHexahedraVerticesVTK
  (std::ostream & out, const VTYPE * hex_vert, const NTYPE numh)
  {
    const NTYPE NUM_VERT_PER_HEXAHEDRON(8);
    using std::endl;

    for (NTYPE i = 0; i < numh; i++) {
      const VTYPE * hex_i_vert = hex_vert + i*NUM_VERT_PER_HEXAHEDRON;

      out << NUM_VERT_PER_HEXAHEDRON;
      out << " " << hex_i_vert[0];
      out << " " << hex_i_vert[1];
      out << " " << hex_i_vert[3];
      out << " " << hex_i_vert[2];
      out << " " << hex_i_vert[4];
      out << " " << hex_i_vert[5];
      out << " " << hex_i_vert[7];
      out << " " << hex_i_vert[6];
      out << endl;
    }
  }


  /// Output scalar data.
  template <typename DATA_TYPE, typename NTYPE> 
  void ijkoutScalarValues
  (std::ostream & out, const DATA_TYPE * data, const NTYPE num)
  {
    using std::endl;

    for (NTYPE i = 0; i < num; i++) 
      {  out << data[i] << endl;  }
  }

  /// Output vtk scalar cell data.
  template <typename DATA_TYPE, typename NTYPE> 
  void ijkoutScalarVTK
  (std::ostream & out, const char * table_name, 
   const char * data_type_string,
   const DATA_TYPE * data, const NTYPE num_cells)
  {
    using std::endl;

    out << "SCALARS " << table_name << " " 
        << data_type_string << " 1" << endl;
    out << "LOOKUP_TABLE default" << endl;

    ijkoutScalarValues(out, data, num_cells);
  }


  /// Output hexahedra.
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1>
  void ijkoutHexahedraVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * hexahedra_vert, const NTYPE1 numh,
   const bool flag_reorder_hex_vertices)
  {
    const int NUM_VERT_PER_HEXAHEDRON(8);
    IJK::PROCEDURE_ERROR error("ijkoutHexahedraVTK");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 available for .vtk files.");
      throw error;
    }

    using std::endl;

    ijkoutVTKheader(out, dataset_name, dim);

    out << "DATASET UNSTRUCTURED_GRID" << endl;
    out << endl;

    out << "POINTS " << numv << " float" << endl;
    ijkoutVertexCoord(out, dim, coord, numv);
    out << endl;

    out << "CELLS " << numh << " " 
        << numh*(1+NUM_VERT_PER_HEXAHEDRON) << endl;

    if (flag_reorder_hex_vertices) {
      ijkoutHexahedraVerticesVTK(out, hexahedra_vert, numh);
    }
    else {
      ijkoutPolygonVertices
        (out, NUM_VERT_PER_HEXAHEDRON, hexahedra_vert, numh);
    }
    out << endl;

    ijkoutHexahedraTypesVTK(out, numh);
  }


  /// Output hexahedra.
  /// - C++ STL vector format for coord[] and hexahedra_vert[].
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE>
  void ijkoutHexahedraVTK
  (std::ostream & out, const char * dataset_name, const int dim,
   const std::vector<CTYPE> & coord, 
   const std::vector<VTYPE> & hexahedra_vert,
   const bool flag_reorder_hex_vertices)
  {
    typedef typename std::vector<CTYPE>::size_type SIZEC_TYPE;
    typedef typename std::vector<VTYPE>::size_type SIZEV_TYPE;

    const SIZEV_TYPE NUM_VERT_PER_HEXAHEDRON(8);
    const SIZEC_TYPE numc = coord.size()/dim;
    const SIZEV_TYPE num_hex = hexahedra_vert.size()/NUM_VERT_PER_HEXAHEDRON;

    ijkoutHexahedraVTK
      (out, dataset_name, dim, vector2pointer(coord), numc,
       vector2pointer(hexahedra_vert), num_hex, flag_reorder_hex_vertices);
  }


  /// Output hexahedra and data for each hex.
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1,
            typename DATA_TYPE>
  void ijkoutHexahedraAndHexDataVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * hexahedra_vert, const NTYPE1 numh,
   const bool flag_reorder_hex_vertices,
   const char * hex_table_name,
   const char * data_type_name,
   const DATA_TYPE * hex_data)
  {
    using namespace std;

    ijkoutHexahedraVTK
      (out, dataset_name, dim, coord, numv, hexahedra_vert, numh,
       flag_reorder_hex_vertices);

    out << endl;

    out << "CELL_DATA " << numh << endl;
    ijkoutScalarVTK
      (out, hex_table_name, data_type_name, hex_data, numh);
  }


  /// Output hexahedra and data for each hex.
  /// - Version where hex_data has type float.
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1>
  void ijkoutHexahedraAndHexDataVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * hexahedra_vert, const NTYPE1 numh,
   const bool flag_reorder_hex_vertices,
   const char * hex_table_name,
   const float * hex_data)
  {
    ijkoutHexahedraAndHexDataVTK
      (out, dataset_name, dim, coord, numv, hexahedra_vert, numh,
       flag_reorder_hex_vertices, hex_table_name, "float", hex_data);
  }


  /// Output hexahedra and data for each hex.
  /// - Version where hex_data is a C++ vector of floats.
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1>
  void ijkoutHexahedraAndHexDataVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * hexahedra_vert, const NTYPE1 numh,
   const bool flag_reorder_hex_vertices,
   const char * hex_table_name,
   const std::vector<float> & hex_data)
  {
    ijkoutHexahedraAndHexDataVTK
      (out, dataset_name, dim, coord, numv, hexahedra_vert, numh,
       flag_reorder_hex_vertices, hex_table_name, vector2pointer(hex_data));
  }


  /// Output hexahedra and data for each hex.
  /// - Version with two arrays of hex_data.
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1,
            typename DATA_A_TYPE, typename DATA_B_TYPE>
  void ijkoutHexahedraAndHexDataVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * hexahedra_vert, const NTYPE1 numh,
   const bool flag_reorder_hex_vertices,
   const char * hex_tableA_name,
   const char * data_typeA_name,
   const DATA_A_TYPE * hex_dataA,
   const char * hex_tableB_name,
   const char * data_typeB_name,
   const DATA_B_TYPE * hex_dataB)
  {
    using namespace std;

    ijkoutHexahedraVTK
      (out, dataset_name, dim, coord, numv, hexahedra_vert, numh,
       flag_reorder_hex_vertices);

    out << "CELL_DATA " << numh << endl;
    ijkoutScalarVTK
      (out, hex_tableA_name, data_typeA_name, hex_dataA, numh);

    out << endl;

    ijkoutScalarVTK
      (out, hex_tableB_name, data_typeB_name, hex_dataB, numh);
  }


  /// Output hexahedra and data for each hex.
  /// - Version with two arrays of hex_data.
  /// - Version where both arrays have type float.
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1>
  void ijkoutHexahedraAndHexDataVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * hexahedra_vert, const NTYPE1 numh,
   const bool flag_reorder_hex_vertices,
   const char * hex_tableA_name,
   const float * hex_dataA,
   const char * hex_tableB_name,
   const float * hex_dataB)
  {
    ijkoutHexahedraAndHexDataVTK
      (out, dataset_name, dim, coord, numv, hexahedra_vert, numh,
       flag_reorder_hex_vertices, 
       hex_tableA_name, "float", hex_dataA, 
       hex_tableB_name, "float", hex_dataB);
  }


  /// Output hexahedra and data for each hex.
  /// - Version with two arrays of hex_data.
  /// - Version where both arrays are C++ vectors of floats.
  /// @param flag_reorder_hex_vertices If true, reorder hex vertices
  ///    in order expected by VTK.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1>
  void ijkoutHexahedraAndHexDataVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * hexahedra_vert, const NTYPE1 numh,
   const bool flag_reorder_hex_vertices,
   const char * hex_tableA_name,
   const std::vector<float> & hex_dataA,
   const char * hex_tableB_name,
   const std::vector<float> & hex_dataB)
  {
    ijkoutHexahedraAndHexDataVTK
      (out, dataset_name, dim, coord, numv, hexahedra_vert, numh,
       flag_reorder_hex_vertices, 
       hex_tableA_name, vector2pointer(hex_dataA),
       hex_tableB_name, vector2pointer(hex_dataB));
  }


  /// Output tetrahedra.
  template <typename CTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1>
  void ijkoutTetrahedraVTK
  (std::ostream & out, const char * dataset_name,
   const int dim, const CTYPE * coord, const NTYPE0 numv,
   const VTYPE * tetrahedra_vert, const NTYPE1 num_tet)
  {
    const int NUM_VERT_PER_TETRAHEDRON(4);
    IJK::PROCEDURE_ERROR error("ijkoutTetrahedraVTK");

    if (dim != 3) {
      error.AddMessage
        ("Programming error.  Only dimension 3 available for .vtk files.");
      throw error;
    }

    using std::endl;

    ijkoutVTKheader(out, dataset_name, dim);

    out << "DATASET UNSTRUCTURED_GRID" << endl;
    out << endl;

    out << "POINTS " << numv << " float" << endl;
    ijkoutVertexCoord(out, dim, coord, numv);
    out << endl;

    out << "CELLS " << num_tet << " " 
        << num_tet*(1+NUM_VERT_PER_TETRAHEDRON) << endl;

    ijkoutPolygonVertices
      (out, NUM_VERT_PER_TETRAHEDRON, tetrahedra_vert, num_tet);
    out << endl;

    ijkoutTetrahedraTypesVTK(out, num_tet);
  }


  /// Output tetrahedra
  /// - C++ STL vector format for coord[] and tetrahedra_vert[].
  template <typename CTYPE, typename VTYPE>
  void ijkoutTetrahedraVTK
  (std::ostream & out, const char * dataset_name, const int dim,
   const std::vector<CTYPE> & coord, 
   const std::vector<VTYPE> & tetrahedra_vert)
  {
    typedef typename std::vector<CTYPE>::size_type SIZEC_TYPE;
    typedef typename std::vector<VTYPE>::size_type SIZEV_TYPE;

    const SIZEV_TYPE NUM_VERT_PER_TETRAHEDRON(4);
    const SIZEC_TYPE numc = coord.size()/dim;
    const SIZEV_TYPE num_tet = tetrahedra_vert.size()/NUM_VERT_PER_TETRAHEDRON;

    ijkoutTetrahedraVTK
      (out, dataset_name, dim, vector2pointer(coord), numc,
       vector2pointer(tetrahedra_vert), num_tet);
  }


  // ******************************************
  // Fig file
  // ******************************************

  // local namespace
  namespace {
    inline int index_other_endpoint(int k)
    {
      int kseg = k/2;
      int j = (k+1)%2;
      return(2*kseg+j);
    }
  }

  // Output Fig polygonal lines.
  template <typename CTYPE, typename VTYPE, typename SCALE_TYPE> 
  void ijkoutFIGpolyline
  (std::ostream & out, const CTYPE * coord, const int numv,
   const VTYPE * seg, const int nums, const SCALE_TYPE scale_coord)
  {
    using std::vector;
    const int cap_style = 1;             // round cap
    IJK::ARRAY<bool> is_edge_processed(nums, false);

    // degree of vertex iv
    IJK::ARRAY<int> degree(numv, 0);

    // next[2*i] = next element (in next[]) around vertex 0 of segment i
    // next[2*i+1] = next element (in next[]) around vertex 1 of segment i
    IJK::ARRAY<int> next(2*nums);

    // first[iv] = first element in next[] incident on vertex iv
    // last[iv] = last element in next[] incident on vertex iv
    IJK::ARRAY<int> first(numv);
    IJK::ARRAY<int> last(numv);

    // set degree[], first[], last[], next[]
    for (int i = 0; i < nums; i++)
      for (int j = 0; j < 2; j++) {
        VTYPE iv = seg[2*i+j];
        if (degree[iv] == 0) {
          first[iv] = 2*i+j;
          last[iv] = 2*i+j;
          next[2*i+j] = 2*i+j;
        }
        else {
          next[2*i+j] = next[last[iv]];
          next[last[iv]] = 2*i+j;
          last[iv] = 2*i+j;
        }
        degree[iv]++;
      }

    vector<VTYPE> vlist;
    for (int i = 0; i < nums; i++) {
      if (is_edge_processed[i]) { continue; };

      vlist.clear();
      VTYPE iv0 = seg[2*i];
      VTYPE iv1 = seg[2*i+1];
      vlist.push_back(iv0);
      vlist.push_back(iv1);
      is_edge_processed[i] = true;

      iv0 = iv1;
      int k = next[2*i+1];
      k = index_other_endpoint(k);
      while (!is_edge_processed[k/2]) {
        iv1 = seg[k];
        vlist.push_back(iv1);
        is_edge_processed[k/2] = true;
        iv0 = iv1;
        k = next[k];
        k = index_other_endpoint(k);
      }

      out << "2 1 0 1 0 7 50 -1 -1 0.000 0 "
          << cap_style << " -1 0 0 " << vlist.size() << std::endl;

      out << "    ";
      for (int j = 0; j < vlist.size(); j++) {
        int iv = vlist[j];
        for (int k = 0; k < 2; k++) {
          int x = int(scale_coord*coord[2*iv+k]);
          out << "  " << x;
        }
      }
      out << std::endl;
    }

  }


  /// Output FIG file header.
  inline void ijkoutFIGheader(std::ostream & out, const bool flag_metric)
  {
    // FIG header
    out << "#FIG 3.2" << std::endl;

    out << "Landscape" << std::endl;
    out << "Center" << std::endl;
    if (flag_metric) { out << "Metric" << std::endl; }
    else { out << "Inches" << std::endl; }
    out << "Letter" << std::endl;
    out << "100.00" << std::endl;
    out << "Single" << std::endl;
    out << "-2" << std::endl;
    out << "1200 2" << std::endl;
  }

  /// \brief Output Fig .fig file.
  /// @param out = output stream
  /// @param dim = dimension.  Must be 2.
  /// @param coord[2*i+k] = k'th coordinate of vertex j  (k < dim)
  /// @param numv = number of vertices 
  /// @param seg[2*j+k] = k'th vertex index of segment j (k < dim)
  /// @param nums = number of line segments
  /// @param flag_polyline: if true, join segments to form FIG polylines
  ///                if false, output segments individually
  /// @param flag_metric: if true, output "Metric" as units
  ///              if false, output "Inches"
  template <typename CTYPE, typename VTYPE, typename SCALE_TYPE> void ijkoutFIG
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * seg, const int nums, const SCALE_TYPE scale_coord,
   const bool flag_polyline = false, const bool flag_metric = false)
  {
    IJK::PROCEDURE_ERROR error("ijkoutFIG");
    const int cap_style = 1;             // round cap

    if (dim != 2)
      throw error("Illegal dimension.  Fig files are only for dimension 2.");

    ijkoutFIGheader(out, flag_metric);

    if (!flag_polyline) {
      // output each line segment separately line segments

      for (int i = 0; i < nums; i++) {
        out << "2 1 0 1 0 7 50 -1 -1 0.000 0 "
            << cap_style << " -1 0 0 2" << std::endl;
        out << "    ";
        for (int j = 0; j < 2; j++) {
          int iv = seg[2*i+j];
          for (int k = 0; k < 2; k++) {
            int x = int(scale_coord*coord[2*iv+k]);
            out << "  " << x;
          }
        }
        out << std::endl;
      }
    }
    else {
      ijkoutFIGpolyline(out, coord, numv, seg, nums, scale_coord);
    }

  }


  /// \brief Output Fig .fig file of polygons.
  /// @param out = output stream
  /// @param dim = dimension.  Must be 2.
  /// @param coord[2*i+k] = k'th coordinate of vertex j  (k < dim)
  /// @param numv = number of vertices 
  /// @param flag_metric: if true, output "Metric" as units
  ///              if false, output "Inches"
  template <typename DTYPE, typename CTYPE, 
            typename NTYPE0, typename NTYPE1, typename NTYPE2,
            typename VTYPE, typename ITYPE,
            typename SCALE_TYPE> void ijkoutPolygonFIG
  (std::ostream & out, const DTYPE dim, 
   const CTYPE * coord, const NTYPE0 numv,
   const NTYPE1 * num_poly_vert, const VTYPE * poly_vert,
   const ITYPE * first_poly_vert, const NTYPE2 num_poly,
   const SCALE_TYPE scale_coord,
   const bool flag_polyline = false, const bool flag_metric = false)
  {
    IJK::PROCEDURE_ERROR error("ijkoutFIG");
    const int cap_style = 1;             // round cap

    if (dim != 2)
      throw error("Illegal dimension.  Fig files are only for dimension 2.");

    ijkoutFIGheader(out, flag_metric);

    for (NTYPE2 ipoly = 0; ipoly < num_poly; ipoly++) {

      out << "2 1 0 1 0 7 50 -1 -1 0.000 0 "
          << cap_style << " -1 0 0 ";

      if (num_poly_vert[ipoly] > 2) {
        out << num_poly_vert[ipoly]+1 << std::endl;

      }
      else {
        out << num_poly_vert[ipoly] << std::endl;
      }

      out << "    ";
      for (NTYPE1 j = 0; j < num_poly_vert[ipoly]; j++) {
        NTYPE1 k = first_poly_vert[ipoly] + j;
        VTYPE iv = poly_vert[k];
        for (int k2 = 0; k2 < 2; k2++) {
          int x = int(scale_coord*coord[2*iv+k2]);
          out << "  " << x;
        }
      }

      if (num_poly_vert[ipoly] > 2) {
        NTYPE1 k = first_poly_vert[ipoly];
        VTYPE iv = poly_vert[k];
        for (int k2 = 0; k2 < 2; k2++) {
          int x = int(scale_coord*coord[2*iv+k2]);
          out << "  " << x;
        }
      }
      out << std::endl;
    }

  }


  // ******************************************
  // Miscellaneous local output routines
  // ******************************************

  namespace {

    /// Output vertex coordinates
    /// @param out Output stream.
    /// @param dim Vertex dimension (number of vertex coordinates.)
    /// @param coord[] Array of vertex coordinates.
    ///                coord[dim*i+k] = k'th coordinate of vertex j  (k < dim).
    /// @param numv Number of vertices.
    template <typename CTYPE> void ijkoutVertexCoord
    (std::ostream & out, const int dim, const CTYPE * coord, const int numv)
    {
      for (int iv = 0; iv < numv; iv++) {
        for (int d = 0; d < dim; d++) {
          out << coord[iv*dim + d];
          if (d < dim-1) { out << " "; }
          else { out << std::endl; };
        }
      }
    }

    /// Output vertex coordinates with color.
    /// @param out Output stream.
    /// @param dim Vertex dimension (number of vertex coordinates.)
    /// @param coord[] Array of vertex coordinates.
    ///                coord[dim*i+k] = k'th coordinate of vertex j  (k < dim).
    /// @param numv Number of vertices.
    template <typename CTYPE, typename COLOR_TYPE> 
    void ijkoutVertexCoordColor
    (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
    {
      for (int iv = 0; iv < numv; iv++) {
        for (int d = 0; d < dim; d++) 
          { out << coord[iv*dim + d] << " "; }

        out << " ";  // Add another space
        
        // Add color
        ijkoutColor(out, iv, front_color, back_color);
        out << std:: endl;
      }
    }

    /// Output polygon vertices
    /// @param out Output stream.
    /// @param numv_per_polygon Number of vertices per polygon.
    /// @pre     All polygons have the same number of vertices.
    /// @param poly_vert[] Array of simplex vertices.
    ///        poly_vert[dim*j+k] = k'th vertex index of polygon j.
    /// @param nump Number of polygons.
    template <typename VTYPE> void ijkoutPolygonVertices
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump)
    {
      for (int is = 0; is < nump; is++) {
        out << numv_per_polygon << " ";
        for (int iv = 0; iv < numv_per_polygon; iv++) {
          out << poly_vert[is*numv_per_polygon + iv];
          if (iv < numv_per_polygon-1) { out << " "; }
          else { out << std::endl; };
        }
      }
    }

    /// Output polytope vertices
    template <typename NTYPE, typename VTYPE, typename ITYPE>
    void ijkoutPolytopeVertices
    (std::ostream & out, 
     const NTYPE * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const int num_poly)
    {
      for (int ipoly = 0; ipoly < num_poly; ipoly++) {
        const ITYPE * pvert = poly_vert + first_poly_vert[ipoly];
        NTYPE num_pvert = num_poly_vert[ipoly];

        out << num_pvert << " ";
        for (int i = 0; i < num_pvert; i++) {
          out << pvert[i];
          if (i+1 < num_pvert) { out << " "; }
          else { out << std::endl; };
        }
      }

    }

    /// Output rgba front and back color.
    template <typename ITYPE, typename COLOR_TYPE>
    void ijkoutColor
    (std::ostream & out, const ITYPE i, 
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
    {
      for (int ic = 0; ic < 3; ic++) 
        { out << front_color[4*i+ic] << " "; }
      out << front_color[4*i+3];

      if (back_color != NULL) {
        out << "  ";
        for (int ic = 0; ic < 3; ic++)
          { out << back_color[4*i+ic] << " "; }
        out << back_color[4*i+3];
      }
    }

    /// Output rgb front and back.
    template <typename ITYPE>
    void ijkoutRGB
    (std::ostream & out, const ITYPE i, 
     const unsigned char * front_rgb, const unsigned char * back_rgb)
    {
      out << int(front_rgb[3*i]);
      out << " " << int(front_rgb[3*i+1]);
      out << " " << int(front_rgb[3*i+2]);

      if (back_rgb != NULL) {
        out << "  ";
        out << int(back_rgb[3*i]);
        out << " " << int(back_rgb[3*i+1]);
        out << " " << int(back_rgb[3*i+2]);
      }
    }

    /// Output rgb.  (No back color.)
    template <typename ITYPE>
    void ijkoutRGB
    (std::ostream & out, const ITYPE i, 
     const unsigned char * front_rgb)
    {
      typedef unsigned char * RGB_PTR_TYPE;

      ijkoutRGB(out, i, front_rgb, RGB_PTR_TYPE(NULL));
    }

    /// Output polygon vertices and polygon color.
    /// @param out Output stream.
    /// @param numv_per_polygon Number of vertices per polygon.
    /// @pre     All polygons have the same number of vertices.
    /// @param poly_vert[] Array of simplex vertices.
    ///        poly_vert[dim*j+k] = k'th vertex index of polygon j.
    /// @param nump Number of polygons.
    template <typename VTYPE, typename COLOR_TYPE> 
    void ijkoutPolygonVerticesColor
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
    {
      for (int is = 0; is < nump; is++) {
        out << numv_per_polygon << " ";
        for (int iv = 0; iv < numv_per_polygon; iv++) {
          out << poly_vert[is*numv_per_polygon + iv] << " ";
        }
        out << " ";  // Add another space
        ijkoutColor(out, is, front_color, back_color);
        out << std::endl;
      }
    }


    /// Output polygon vertices and polygon color.
    /// Single color per polygon.
    /// @param out Output stream.
    /// @param numv_per_polygon Number of vertices per polygon.
    /// @pre     All polygons have the same number of vertices.
    /// @param poly_vert[] Array of simplex vertices.
    ///        poly_vert[dim*j+k] = k'th vertex index of polygon j.
    /// @param nump Number of polygons.
    template <typename VTYPE, typename COLOR_TYPE> 
    void ijkoutPolygonVerticesColor
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump,
     const COLOR_TYPE * poly_color)
    {
      for (int is = 0; is < nump; is++) {
        out << numv_per_polygon << " ";
        for (int iv = 0; iv < numv_per_polygon; iv++) {
          out << poly_vert[is*numv_per_polygon + iv] << " ";
        }
        out << " ";  // Add another space
        ijkoutColor(out, is, poly_color);
        out << std::endl;
      }
    }

    /// Output polygon vertices and polygon rgb.
    /// Single rgb per polygon.
    /// @param out Output stream.
    /// @param numv_per_polygon Number of vertices per polygon.
    /// @pre     All polygons have the same number of vertices.
    /// @param poly_vert[] Array of simplex vertices.
    ///        poly_vert[dim*j+k] = k'th vertex index of polygon j.
    /// @param nump Number of polygons.
    template <typename VTYPE> 
    void ijkoutPolygonVerticesRGB
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump,
     const unsigned char * poly_rgb)
    {
      for (int is = 0; is < nump; is++) {
        out << numv_per_polygon << " ";
        for (int iv = 0; iv < numv_per_polygon; iv++) {
          out << poly_vert[is*numv_per_polygon + iv] << " ";
        }
        out << " ";  // Add another space

        ijkoutRGB(out, is, poly_rgb);
        out << std::endl;
      }
    }


    /// Output polygon vertices and polygon rgb.
    /// - Version where polygons may have different numbers of vertices.
    /// - Single rgb per polygon.
    /// @param out Output stream.
    /// @param numv_per_polygon Number of vertices per polygon.
    /// @pre     All polygons have the same number of vertices.
    /// @param poly_vert[] Array of simplex vertices.
    ///        poly_vert[dim*j+k] = k'th vertex index of polygon j.
    /// @param nump Number of polygons.
    template <typename NTYPE0, typename NTYPE1,
              typename VTYPE, typename ITYPE>
    void ijkoutPolygonVerticesRGBII
    (std::ostream & out, 
     const NTYPE0 * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const NTYPE1 num_poly,
     const unsigned char * poly_rgb)
    {
      for (NTYPE1 ipoly = 0; ipoly < num_poly; ipoly++) {
        const NTYPE0 numv = num_poly_vert[ipoly];
        out << numv << " ";
        for (int iv = 0; iv < numv; iv++) {
          out << poly_vert[ipoly*numv + iv] << " ";
        }
        out << " ";  // Add another space

        ijkoutRGB(out, ipoly, poly_rgb);
        out << std::endl;
      }
    }


    /// Output polytope vertices and polytope color.
    /// Different polytopes may have different numbers of vertices.
    /// @param out Output stream.
    template <typename NTYPE, typename VTYPE, typename ITYPE,
              typename COLOR_TYPE>
    void ijkoutPolytopeVerticesColor
    (std::ostream & out, 
     const NTYPE * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const int num_poly,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
    {
      for (int ipoly = 0; ipoly < num_poly; ipoly++) {
        const ITYPE * pvert = poly_vert + first_poly_vert[ipoly];
        NTYPE num_pvert = num_poly_vert[ipoly];

        out << num_pvert << " ";
        for (int i = 0; i < num_pvert; i++) {
          out << pvert[i] << " ";
        }

        out << " ";  // Add another space
        
        // Add color
        ijkoutColor(out, ipoly, front_color, back_color);
        out << std:: endl;
      }
    }

    /// Output quad vertices.
    /// @param out Output stream.
    /// @param quad_vert[] Array of quad vertices.
    ///        quad_vert[4*j+k] = k'th vertex index of polygon j.
    /// @param numq Number of quadrilaterals.
    ///  @param flag_reorder_vertices = if true, change vertex order to be
    ///                         counter-clockwise around quad.
    template <typename VTYPE> void ijkoutQuadVertices
    (std::ostream & out, const VTYPE * quad_vert, const int numq,
     const bool flag_reorder_vertices)
    {
      const int NUMV_PER_QUAD = 4;

      if (flag_reorder_vertices) {
        for (int iq = 0; iq < numq; iq++) {
          out << NUMV_PER_QUAD << " ";
          const VTYPE * v = quad_vert+iq*NUMV_PER_QUAD;
          out << v[0] << " " << v[1] << " ";

          // Note change in order between v[2] and v[3]
          out << v[3] << " " << v[2] << std::endl;
        }
      }
      else {
        ijkoutPolygonVertices(out, NUMV_PER_QUAD, quad_vert, numq);
      }

    }

    /// Output vector endpoint coordinates.
    /// @param out Output stream.
    template <typename CTYPE0, typename CTYPE1, typename CTYPE2> 
    void ijkoutVectorEndpoints
    (std::ostream & out, const int dim, 
     const CTYPE0 * point_coord, const CTYPE1 * dir_coord, 
     const int num_vectors, const CTYPE2 scale)
    {
      for (int iv = 0; iv < num_vectors; iv++) {

        // Output point iv.
        for (int d = 0; d < dim; d++) {
          out << point_coord[iv*dim + d];
          if (d < dim-1) { out << " "; }
          else { out << std::endl; };
        }

        // Output coordinates of (point iv)+(direction iv)
        for (int d = 0; d < dim; d++) {
          int j = iv*dim + d;
          out << (point_coord[j] + dir_coord[j]*scale);
          if (d < dim-1) { out << " "; }
          else { out << std::endl; };
        }
      }
    }

  }

}

#endif
