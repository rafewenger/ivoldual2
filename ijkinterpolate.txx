/// \file ijkinterpolate.txx
/// ijk templates for linear and multilinear interpolation
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009-2017 Rephael Wenger

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

#ifndef _IJKINTERPOLATE_
#define _IJKINTERPOLATE_

#include "ijkcoord.txx"

namespace IJK {

  // **********************************************************
  //! @name Linear and multilinear interpolation functions.
  // **********************************************************

  ///@{

  /// Compute coordinates using linear interpolation.
  template <typename DTYPE, typename WTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void linear_interpolate_coord
  (const DTYPE dimension, const WTYPE w0, const CTYPE0 coord0[],
   const CTYPE1 coord1[], CTYPE2 coord2[])
  {
    const WTYPE w1 = 1-w0;

    for (int d = 0; d < dimension; d++)
      { coord2[d] = CTYPE2(w0*coord0[d] + w1*coord1[d]); }
  }

  /// Compute coordinates using linear interpolation.
  /// @param PTYPE = Precision type to be used in calculations. 
  ///                Should be float or double.
  template <typename PTYPE, typename DTYPE, typename STYPE, typename STYPE2,
            typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void linear_interpolate_coord
  (const DTYPE dimension, const STYPE s0, const CTYPE0 coord0[],
   const STYPE s1, const CTYPE1 coord1[],
   const STYPE2 isovalue, CTYPE2 coord2[])
  {
    PTYPE w0;
    STYPE s_diff = s1 - s0;
    const PTYPE EPSILON = 0.00001;
    if (s_diff > EPSILON || s_diff < -EPSILON) { 
      w0 = (s1 - isovalue) / s_diff;
    }
    else {
      // arbitrarily set weights to 0.5
      w0 = 0.5;
    };

    linear_interpolate_coord
      (dimension, w0, coord0, coord1, coord2);
  }

  /// Compute coordinates using linear interpolation.
  template <typename DTYPE, typename STYPE, typename VTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void linear_interpolate_coord
  (const DTYPE dimension, const STYPE s0, const CTYPE0 coord0[],
   const STYPE s1, const CTYPE1 coord1[],
   const VTYPE isovalue, CTYPE2 coord2[])
  {
    linear_interpolate_coord<double>
      (dimension, s0, coord0, s1, coord1, isovalue, coord2);
  }

  /// Multilinear interpolate scalar values at unit cube vertices.
  /// @tparam ITYPE = Type of interpolated scalar value.
  /// @tparam WTYPE = Weight type.
  template <typename ITYPE, typename WTYPE,
            typename DTYPE, typename CTYPE, typename NTYPE, typename STYPE>
  ITYPE multilinear_interpolate_scalar
  (const DTYPE dimension, const CTYPE coord[],
   const NTYPE num_cube_vertices, const STYPE cube_scalar[])
  {
    ITYPE s = 0;
    for (NTYPE i = 0; i < num_cube_vertices; i++) {
      WTYPE weight = 1.0;
      NTYPE mask = 1;
      for (DTYPE d = 0; d < dimension; d++) {
        if ((i & mask) == 0) { weight = weight * (1-coord[d]); }
        else { weight = weight*coord[d]; };

        mask = (mask << 1L);
      }
      s += weight * cube_scalar[i];
    }
    return(s);
  }

  /// Multilinear interpolate scalar values at unit cube vertices.
  template <typename DTYPE, typename CTYPE, typename NTYPE, typename STYPE>
  double multilinear_interpolate_scalar
  (const DTYPE dimension, const CTYPE coord[],
   const NTYPE num_cube_vertices, const STYPE cube_scalar[])
  {
    double s = multilinear_interpolate_scalar<double, double>
      (dimension, coord, num_cube_vertices, cube_scalar);

    return(s);
  }

  namespace {

    template <typename STYPE, typename VTYPE>
    inline bool get_vertex_sign
    (const STYPE s, const VTYPE isovalue)
    {
      if (s < isovalue) { return(false); }
      else { return(true); };
    }

  }


  /// Compute location of isosurface vertex on line segment [A,B].
  /// Use multilinear interpolation on unit cube vertices 
  ///   to calculate scalar values.
  /// @tparam ITYPE = Type of interpolated scalar value.
  /// @tparam WTYPE = Weight type.
  template <typename ITYPE, typename WTYPE, typename DTYPE, 
            typename CTYPE, typename VTYPE, typename NTYPE, typename STYPE>
  void multilinear_interpolate_coord
  (const DTYPE dimension, const CTYPE coordA[], const CTYPE coordB[],
   const VTYPE isovalue, const NTYPE num_cube_vertices,
   const STYPE cube_scalar[], const NTYPE num_iter, CTYPE coordC[])
  {

    ITYPE s0 = multilinear_interpolate_scalar<ITYPE,WTYPE>
      (dimension, coordA, num_cube_vertices, cube_scalar);
    ITYPE s1 = multilinear_interpolate_scalar<ITYPE,WTYPE>
      (dimension, coordB, num_cube_vertices, cube_scalar);

    bool sign0 = get_vertex_sign(s0, isovalue);
    bool sign1 = get_vertex_sign(s1, isovalue);

    if (sign0 == sign1) {
      linear_interpolate_coord(dimension, 0.5, coordA, coordB, coordC);
      return;
    }

    if (s0 == isovalue) {
      IJK::copy_coord(dimension, coordA, coordC);
      return;
    }

    if (s1 == isovalue) {
      IJK::copy_coord(dimension, coordB, coordC);
      return;
    }

    WTYPE t0 = 0.0;
    WTYPE t1 = 1.0;
    WTYPE t2;
    ITYPE s2;
    for (NTYPE i = 0; i  < num_iter; i++) {
      t2 = (t0+t1)/2.0;
      linear_interpolate_coord(dimension, 1.0-t2, coordA, coordB, coordC);
      s2 = multilinear_interpolate_scalar<ITYPE,WTYPE>
        (dimension, coordC, num_cube_vertices, cube_scalar);
      bool sign2 = get_vertex_sign(s2, isovalue);

      if (sign0 == sign2) { t0 = t2; }
      else { t1 = t2; }
    }

  }


  /// Compute location of isosurface vertex on line segment [A,B].
  /// Use multilinear interpolation on unit cube vertices 
  ///   to calculate scalar values.
  template <typename DTYPE, 
            typename CTYPE, typename VTYPE, typename NTYPE, typename STYPE>
  void multilinear_interpolate_coord
  (const DTYPE dimension, const CTYPE coordA[], const CTYPE coordB[],
   const VTYPE isovalue, const NTYPE num_cube_vertices,
   const STYPE cube_scalar[], const NTYPE num_iter, CTYPE coordC[])
  {
    multilinear_interpolate_coord<double,double>
      (dimension, coordA, coordB, isovalue, num_cube_vertices,
       cube_scalar, num_iter, coordC);
  }

  ///@}


  // **********************************************************
  //! @name Compute centroid.
  // **********************************************************

  ///@{

  /// Compute centroid of a list of points.
  /// @param point_list[] List of points.  Indices into array coord[].
  /// @param coord[] List of point coordinates.
  ///   i'th point has coordinates starting at coord[point_list[i]*dimension].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, typename PTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE>
  void compute_coord_centroid
  (const DTYPE dimension, const PTYPE point_list[], const NTYPE num_points,
   const CTYPE0 coord[], CTYPE1 centroid_coord[])
  {
    IJK::set_coord(dimension, 0, centroid_coord);
    
    for (NTYPE i = 0; i < num_points; i++) {
      PTYPE ipoint = point_list[i];
      IJK::add_coord
        (dimension, coord+ipoint*dimension, centroid_coord, centroid_coord);
    }

    if (num_points > 0) {
      IJK::divide_coord(dimension, num_points, centroid_coord, centroid_coord);
    }
  }


  /// Compute centroid of a list of points.
  /// Version using C++ vector of coordinates.
  template <typename DTYPE, typename PTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE>
  void compute_coord_centroid
  (const DTYPE dimension, const PTYPE point_list[], const NTYPE num_points,
   const std::vector<CTYPE0> & coord, CTYPE1 centroid_coord[])
  {
    compute_coord_centroid
      (dimension, point_list, num_points, IJK::vector2pointer(coord),
       centroid_coord);
  }


  /// Compute centroid of quad.
  template <typename DTYPE, typename ISOV_INDEX_TYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1>
  void compute_quad_centroid
  (const DTYPE dimension, const ISOV_INDEX_TYPE quad_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 centroid_coord[])
  {
    const int NUM_VERT_PER_QUAD(4);

    compute_coord_centroid
      (dimension, quad_vert, NUM_VERT_PER_QUAD, coord_array, centroid_coord);
  }


  /// Compute centroid of a pentagon.
  template <typename DTYPE, typename ISOV_INDEX_TYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1>
  void compute_pentagon_centroid
  (const DTYPE dimension, const ISOV_INDEX_TYPE pentagon_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 centroid_coord[])
  {
    const int NUM_VERT_PER_PENTAGON(5);

    compute_coord_centroid
      (dimension, pentagon_vert, NUM_VERT_PER_PENTAGON, coord_array, 
       centroid_coord);
  }

  /// Compute centroid of hexahedron.
  template <typename DTYPE, typename ISOV_INDEX_TYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1>
  void compute_hexahedron_centroid
  (const DTYPE dimension, const ISOV_INDEX_TYPE hexahedron_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 centroid_coord[])
  {
    const int NUM_VERT_PER_HEXAHEDRON(8);

    compute_coord_centroid
      (dimension, hexahedron_vert, NUM_VERT_PER_HEXAHEDRON, coord_array, 
       centroid_coord);
  }

  ///@}


  // **********************************************************
  //! @name Compute triangle incenter.
  // **********************************************************

  ///@{

  /// Compute incenter of triangle formed by 3 points.
  /// @param triangle_vert[] Triangle vertices.
  /// @param coord[] List of point coordinates.
  ///   i'th point has coordinates starting at coord[point_list[i]*dimension].
  /// @param[out] incenter_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, typename PTYPE, typename CTYPE0, typename CTYPE1>
  void compute_triangle_incenter
  (const DTYPE dimension, const PTYPE triangle_vert[], const CTYPE0 coord[], 
   CTYPE1 incenter_coord[])
  {
    const PTYPE v0 = triangle_vert[0];
    const PTYPE v1 = triangle_vert[1];
    const PTYPE v2 = triangle_vert[2];
    CTYPE0 length01, length12, length02;

    compute_distance
      (dimension, coord+v0*dimension, coord+v1*dimension, length01);
    compute_distance
      (dimension, coord+v1*dimension, coord+v2*dimension, length12);
    compute_distance
      (dimension, coord+v0*dimension, coord+v2*dimension, length02);

    const CTYPE0 perimeter = length01 + length12 + length02;
    const CTYPE0 w2 = length01/perimeter;
    const CTYPE0 w0 = length12/perimeter;
    const CTYPE0 w1 = length02/perimeter;

    compute_weighted_average_3coord
      (dimension, w0, coord+v0*dimension, w1, coord+v1*dimension,
       w2, coord+v2*dimension, incenter_coord);
  }


  /// Compute incenter of triangle formed by 3 points.
  /// Version using C++ vector of coordinates.
  template <typename DTYPE, typename PTYPE, typename CTYPE0, typename CTYPE1>
  void compute_triangle_incenter
  (const DTYPE dimension, const PTYPE triangle_vert[],
   const std::vector<CTYPE0> & coord, CTYPE1 incenter_coord[])
  {
    compute_triangle_incenter
      (dimension, triangle_vert, IJK::vector2pointer(coord), incenter_coord);
  }

  ///@}


  // **********************************************************
  /// @name Intersect line segment and hyperplane.
  // **********************************************************

  ///@{

  /// Compute intersection of line segment and hyperplane
  ///   using linear interpolation.
  /// @pre Line segment intersects hyperplane.
  /// @param coord0[] Coordinates of endpoint 0 of line segment.
  /// @param coord1[] Coordinates of endpoint 1 of line segment.
  /// @param hplane_orth_dir Direction orthogonal to hyperplane.
  /// @param hplane_coord Hyperplane hplane_orth_dir coordinate.
  /// @param intersection_coord[] Intersection coordinates.
  /// @param EPSILON If |coord0[hplane_orth_dir]-coord1[hplane_orth_dir]|
  ///    is less than epsilon, choose line segment midpoint.
  /// @pre EPSILON >= 0.
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, typename CTYPE3,
            typename DIR_TYPE>
  void compute_line_segment_hyperplane_intersection
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[], 
   const DIR_TYPE hplane_orth_dir, const CTYPE2 hplane_coord,
   CTYPE3 intersection_coord[],
   const double EPSILON = 0.0001)
  {
    const CTYPE0 c0 = coord0[hplane_orth_dir];
    const CTYPE1 c1 = coord1[hplane_orth_dir];
    CTYPE3 w0;

    CTYPE3 cdiff = c0 - c1;
     if (cdiff > EPSILON || cdiff < -EPSILON) {
      w0 = (hplane_coord - c1)/(cdiff);
    }
    else {
      w0 = 0.5;
    }

    if (w0 < 0) { w0 = 0; }
    else if (w0 > 1) { w0 = 1; }

    IJK::linear_interpolate_coord
      (dimension, w0, coord0, coord1, intersection_coord);

    intersection_coord[hplane_orth_dir] = hplane_coord;
  }

  /// Compute intersection of line segment and hyperplane
  ///   using linear interpolation.
  /// Version using array of coordinates.
  /// @param coord_array[] Array of point coordinates.
  /// @param ipoint0 Index of first point in coord_array[].
  /// @param ipoint1 Index of second point in coord_array[].
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename ITYPE0, typename ITYPE1, typename DIR_TYPE>
  void compute_line_segment_hyperplane_intersection
  (const DTYPE dimension, const CTYPE0 coord_array[], 
   const ITYPE0 ipoint0, const ITYPE1 ipoint1,
   const DIR_TYPE hplane_orth_dir, const CTYPE1 hplane_coord,
   CTYPE2 intersection_coord[],
   const double EPSILON = 0.0001)
  {
    const CTYPE0 * coord0 = coord_array + dimension*ipoint0;
    const CTYPE0 * coord1 = coord_array + dimension*ipoint1;

    compute_line_segment_hyperplane_intersection
      (dimension, coord0, coord1, hplane_orth_dir, hplane_coord,
       intersection_coord, EPSILON);
  }

  /// Compute intersection of line segment and hyperplane
  ///   using linear interpolation.
  /// Version using C++ vector of coordinates.
  /// @param coord_array[] Array of point coordinates.
  /// @param ipoint0 Index of first point in coord_array[].
  /// @param ipoint1 Index of second point in coord_array[].
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename ITYPE0, typename ITYPE1, typename DIR_TYPE>
  void compute_line_segment_hyperplane_intersection
  (const DTYPE dimension, const std::vector<CTYPE0> coord_array, 
   const ITYPE0 ipoint0, const ITYPE1 ipoint1,
   const DIR_TYPE hplane_orth_dir, const CTYPE1 hplane_coord,
   CTYPE2 intersection_coord[],
   const double EPSILON = 0.0001)
  {
    compute_line_segment_hyperplane_intersection
      (dimension, IJK::vector2pointer(coord_array),
       ipoint0, ipoint1, hplane_orth_dir, hplane_coord,
       intersection_coord, EPSILON);
  }

  ///@}


  // *********************************************************************
  /// @name Calculation of saddle points of multilinear interpolants.
  // *********************************************************************

  ///@{

  /// Compute parity of integer i.  Parity is 0 or 1.
  template <typename DTYPE, typename ITYPE>
  ITYPE compute_parity(const DTYPE dimension, const ITYPE i)
  {
    ITYPE i2 = i;
    ITYPE parity = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      parity = (parity+i2)%2;
      i2 = (i2 >> 1);
    }
    return(parity);
  }

  /// Convert parity to sign. sign = (-1)^parity.
  template <typename ITYPE>
  ITYPE convert_parity_to_sign(const ITYPE parity)
  {
    if (parity == 0) { return(1); }
    else { return(-1); }
  }

  /// Compute coefficient of highest degree term in multilinear representation.
  template <typename CTYPE, typename DTYPE, typename STYPE, typename ITYPE>
  CTYPE compute_multilinear_a
  (const DTYPE dimension, const STYPE * scalar, const ITYPE * parity_sign)
  {
    // number of cube vertices
    const long numv = (1L << dimension);

    CTYPE a = 0;
    for (long i = 0; i < numv; i++) {
      ITYPE sign = convert_parity_to_sign(dimension%2); 
      sign = sign * parity_sign[i];
      a += scalar[i] * sign;
    }

    return(a);
  }

  /// Compute coefficient of second highest degree terms in multilinear representation.
  /// @param j = Index of variable missing in second highest degree term.
  template <typename CTYPE, typename DTYPE, typename STYPE, typename ITYPE>
  CTYPE compute_multilinear_b
  (const DTYPE dimension, const STYPE * scalar, const ITYPE * parity_sign,
   const DTYPE j)
  {
    // number of cube vertices
    const long numv = (1L << dimension);

    CTYPE b = 0;
    long mask = (1L << j);
    for (long i = 0; i < numv; i++) {

      if ((i & mask) == 0) {
        ITYPE sign = convert_parity_to_sign((dimension+1)%2); 
        sign = sign * parity_sign[i];
        b += scalar[i] * sign;
      }
    }

    return(b);
  }

  /// Compute coefficient of third highest degree terms in multilinear representation.
  /// @param j1 = Index of first variable missing in second highest degree term.
  /// @param j2 = Index of second variable missing in second highest degree term.
  template <typename CTYPE, typename DTYPE, typename STYPE, typename ITYPE>
  CTYPE compute_multilinear_c
  (const DTYPE dimension, const STYPE * scalar, const ITYPE * parity_sign,
   const DTYPE j1, const DTYPE j2)
  {
    // number of cube vertices
    const long numv = (1L << dimension);

    CTYPE c = 0;
    long mask = ((1L << j1) | (1L << j2));
    for (long i = 0; i < numv; i++) {

      if ((i & mask) == 0) {
        ITYPE sign = convert_parity_to_sign(dimension%2); 
        sign = sign * parity_sign[i];
        c += scalar[i] * sign;
      }
    }

    return(c);
  }

  /// Compute saddle points of 3D unit cube.
  /// Saddle points are trilinear_mipdoint[] +/- saddle_offset[].
  template <typename STYPE, typename ITYPE, typename ATYPE, typename CTYPE,
            typename NTYPE, typename ETYPE>
  void compute_saddle_3D
  (const STYPE * scalar, const ITYPE * parity_sign,
   ATYPE & a, CTYPE trilinear_midpoint[3],
   CTYPE saddle_offset[3], NTYPE & num_saddle,
   const ETYPE EPSILON)
  {
    const int dimension = 3;
    ATYPE b[3];
    ATYPE c[3];
    CTYPE q[3];

    // initialize
    num_saddle = 0;
    a = 0;
    for (int d = 0; d < dimension; d++) {
      trilinear_midpoint[d] = 0;
      saddle_offset[d] = 0;
    }

    a = compute_multilinear_a<STYPE>(dimension, scalar, parity_sign);

    if (-EPSILON < a && a < EPSILON) { 
      a = 0;
      return; 
    }

    for (int i = 0; i < 3; i++) {
      b[i] = compute_multilinear_b<STYPE>
        (dimension, scalar, parity_sign, i);

      int j1 = (i+1)%dimension;
      int j2 = (i+2)%dimension;
      c[i] = compute_multilinear_c<STYPE>
        (dimension, scalar, parity_sign, j1, j2);
    }

    for (int d = 0; d < 3; d++) 
      { trilinear_midpoint[d] = -b[d]/a; }

    for (int i = 0; i < 3; i++) {
      int j1 = (i+1)%dimension;
      int j2 = (i+2)%dimension;

      q[i] = b[j1]*b[j2] - a * c[i];

      // Round q[i] to zero if it is very close to 0.
      if (-EPSILON < q[i] && q[i] < EPSILON)
        { q[i] = 0; };
    }

    CTYPE t = q[0]*q[1]*q[2];

    if (q[0] == 0 && q[1] == 0 && q[2] == 0) {
      // One saddle.
      num_saddle = 1;
      return;
    }

    if (t <= 0) { 
      // Zero saddles.
      num_saddle = 0;
      return; 
    }

    CTYPE r = 0;
    r = sqrt(t);
    for (int i = 0; i < 3; i++) { saddle_offset[i] = r / (a*q[i]); }
    num_saddle = 2;
  }

  ///@}

}

#endif
