/// \file ijklist.txx
/// ijk templates for handling lists.
/// - Version 0.2.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010-2017 Rephael Wenger

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

#ifndef _IJKLIST_
#define _IJKLIST_

#include "ijk.txx"

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>

namespace IJK {

  // **************************************************
  // REMAP LIST
  // **************************************************

  /// Find "first" (identifying) element of the set containing x.
  /// Modifies set with path compression.
  /// @pre T is an integer type.
  template <typename T1, typename T2>
  T2 find_set(const T1 x, T2 * set_ident)
  {
    // Find first element
    T1 y = set_ident[x];
    while (y != set_ident[y])
      { y = set_ident[y]; }

    // apply path compression
    T1 z = set_ident[x];
    while (z != set_ident[z]) {
      T1 z_next = set_ident[z];
      set_ident[z] = y;
      z = z_next;
    }

    return(y);
  }

  /// Find "first" (identifying) element of the set containing x.
  /// - Modifies set with path compression.
  /// @pre T is an integer type.
  /// @pre set_ident.size() > 0.
  /// - Version using C++ STL vector for array set_ident[].
  template <typename T1, typename T2>
  T2 find_set(const T1 x, std::vector<T2> & set_ident)
  {
    return(find_set(x, &(set_ident.front())));
  }

  /// Remap members of list[] based on pairs in remap_pair[]
  /// @pre T1, T2 and T3 must be an integer type, i.e. short, int, long.
  /// @pre list[] is a list of non-negative integers.
  /// @pre max_element is a non-negative integer greater than or equal to
  ///      any number in list[]
  template <typename T1, typename T2, typename T3, typename T4,
            typename NTYPE>
  void remap_list(const std::vector< std::pair< T1, T2 > > & remap_pair,
                  const T4 max_element, T3 * list, const NTYPE list_length)
  {
    IJK::PROCEDURE_ERROR error("remap_list");

    if (list_length <= 0 || remap_pair.size() == 0) { return; }

    if (max_element < 0) {
      error.AddMessage
        ("Programming error.  max_element must be non-negative.");
      error.AddMessage("  max_element = ", max_element, ".");
      throw error;
    };

    std::vector<T3> remap_set(max_element+1);
    for (T3 i = 0; i <= max_element; i++) 
      { remap_set[i] = i; }

    for (T3 i = 0; i < remap_pair.size(); i++) {
      T1 k1 = remap_pair[i].first;
      T2 k2 = remap_pair[i].second;

      if (k1 <= max_element && k2 <= max_element) {
        remap_set[k1] = find_set(k2, remap_set);
      }
    }

    for (T3 i = 0; i < list_length; i++) 
      { list[i] = find_set(list[i], remap_set); }
  }

  /// Remap members of list[] based on pairs in remap_pair[]
  /// @pre T1, T2 and T3 must be an integer type, i.e. short, int, long.
  /// @pre list[] is a list of non-negative integers.
  template <typename T1, typename T2, typename T3, typename NTYPE>
  void remap_list(const std::vector< std::pair< T1, T2 > > & remap_pair,
                  T3 * list, const NTYPE list_length)
  {
    IJK::PROCEDURE_ERROR error("remap_list");

    if (list_length <= 0 || remap_pair.size() == 0) { return; }

    T3 * max_ptr = std::max_element(list, list+list_length);

    if (*max_ptr < 0) {
      error.AddMessage
        ("Programming error.  Only non-negative numbers allowed in list[].");
      error.AddMessage("  List contains ", *max_ptr, ".");
      throw error;
    };

    remap_list(remap_pair, *max_ptr, list, list_length);
  }

  /// Remap members of list[] based on pairs in remap_pair[]
  /// @pre T1, T2 and T3 must be an integer type, i.e. short, int, long.
  /// @pre list[] is a list of non-negative integers.
  template <typename T1, typename T2, typename T3>
  void remap_list(const std::vector< std::pair< T1, T2 > > & remap_pair,
                  std::vector<T3> & list)
  {
    if (list.size() == 0) { return; }

    remap_list(remap_pair, &(list.front()), list.size());
  }

  // **************************************************
  // GET DISTINCT LIST ELEMENTS
  // **************************************************

  /// Get distinct elements of list[] and store in distinct_list[].
  /// Retains the order of the distinct elements given in list[].
  /// O(n^2) appropriate only for small lists.
  /// @param distinct_list[] is preallocated with length at least list_length.
  template <typename T1, typename T2, typename NTYPE1, typename NTYPE2>
  void get_distinct(const T1 * list, const NTYPE1 list_length,
                    T2 * distinct_list, NTYPE2 & num_distinct)
  {
    if (list_length == 0) { 
      num_distinct = 0;
      return; 
    };

    distinct_list[0] = list[0];
    num_distinct = 1;

    for (NTYPE1 i = 1; i < list_length; i++) {
      const T1 * previous_ptr = std::find(list, list+i, list[i]);

      if (previous_ptr == list+i) {
        // new element
        distinct_list[num_distinct] = list[i];
        num_distinct++;
      }
    }

  }


  // **************************************************
  // QUERY LIST
  // **************************************************


  /// Return true if list contains element el.
  /// @param[out] iloc Location (index) of element el in list.
  /// iloc is undefined if list does not contain element el.
  template <typename T1, typename T2, typename NTYPE1, typename NTYPE2>
  bool does_list_contain
  (const T1 * list, const NTYPE1 list_length, const T2 el, NTYPE2 & iloc)
  {
    for (iloc = 0; iloc < list_length; iloc++) {
      if (list[iloc] == el) { return(true); }
    }

    return(false);
  }


  /// Return true if list contains element el.
  /// - Version which does not return location of iloc.
  template <typename T1, typename T2, typename NTYPE>
  bool does_list_contain
  (const T1 * list, const NTYPE list_length, const T2 el)
  {
    NTYPE iloc;
    return(does_list_contain(list, list_length, el, iloc));
  }


  /// Return true if list contains element el.
  /// - Version which does not return location of iloc.
  /// - Version using C++ STL vector for array list[].
  template <typename T1, typename T2>
  bool does_list_contain
  (const std::vector<T1> & list, const T2 el)
  {
    return(does_list_contain(vector2pointer(list), list.size(), el));
  }


  /// Return true if elements of listA equal elements of listB
  ///   and are in the same order in each list.
  /// @pre Lists listA and listB have the same length.
  template <typename TA, typename TB, typename NTYPE>
  bool are_lists_equal
  (const TA * listA, const TB * listB, const NTYPE list_length)
  {
    for (NTYPE i = 0; i < list_length; i++) {
      if (listA[i] != listB[i]) { return(false); }
    }

    return(true);
  }


  /// Return true if listA and listB contain a common element.
  /// @param[out] ilocA Location (index) of common element in listA.
  /// ilocA is undefined if listA does not contain any element of listB.
  /// @param[out] ilocB Location (index) of common element in listB.
  /// ilocB is undefined if listA does not contain any element of listB.
  template <typename TA, typename TB, typename NTYPEA, typename NTYPEB,
            typename ITYPEA, typename ITYPEB>
  bool do_lists_intersect
  (const TA * listA, const NTYPEA listA_length,
   const TB * listB, const NTYPEB listB_length,
   ITYPEA & ilocA, ITYPEB & ilocB)
  {
    for (ilocB = 0; ilocB < listB_length; ilocB++) {
      if (does_list_contain(listA, listA_length, listB[ilocB], ilocA)) 
        { return(true); }
    }

    return(false);
  }


  /// Return true if listA and listB contain a common element.
  /// - Version which does not return location of ilocA or ilocB.
  template <typename TA, typename TB, typename NTYPEA, typename NTYPEB>
  bool do_lists_intersect
  (const TA * listA, const NTYPEA listA_length,
   const TB * listB, const NTYPEB listB_length)
  {
    NTYPEA ilocA;
    NTYPEB ilocB;

    return(do_lists_intersect
           (listA, listA_length, listB, listB_length, ilocA, ilocB));
  }

  /// Return true if listA and listB contain a common element.
  /// - Version which does not return location of ilocA or ilocB.
  /// - Version where listB is a C++ vector.
  template <typename TA, typename TB, typename NTYPEA>
  bool do_lists_intersect
  (const TA * listA, const NTYPEA listA_length, const std::vector<TB> & listB)
  {
    return(do_lists_intersect
           (listA, listA_length, vector2pointer(listB), listB.size()));
  }


  /// Return true if list[i] = list[i+1] for some i
  ///   or list[0] = list[list_length-1].
  template <typename T1, typename NTYPE>
  bool does_circular_list_contain_adjacent_duplicate
  (const T1 * list, const NTYPE list_length)
  {
    if (list_length == 0) { return(false); }

    if (list[0] == list[list_length-1]) { return(true); }

    for (NTYPE i = 0; i+1 < list_length; i++) {
      if (list[i] == list[i+1]) { return(true); }
    }

    return(false);
  }


  /// Return next element in list after el.
  /// @pre List contains el.
  template <typename T1, typename T2, typename NTYPE>
  T1 next_element
  (const T1 * list, const NTYPE list_length, const T2 el)
  {
    for (NTYPE i = 0; i < list_length; i++) {
      if (list[i] == el) { 
        NTYPE i2 = (i+1)%list_length;
        return(list[i2]); 
      }
    }

    // Just in case el is not in list.
    return(el);
  }

  /// Return previous element in list before el.
  /// @pre List contains el.
  template <typename T1, typename T2, typename NTYPE>
  T1 prev_element
  (const T1 * list, const NTYPE list_length, const T2 el)
  {
    for (NTYPE i = 0; i < list_length; i++) {
      if (list[i] == el) { 
        NTYPE i2 = (i+list_length-1)%list_length;
        return(list[i2]); 
      }
    }

    // Just in case el is not in list.
    return(el);
  }


  // **************************************************
  // ADD LIST
  // **************************************************

  /// Add list to a C++ STL vector.
  template <typename T1, typename T2, typename NTYPE>
  void add_list(const T1 * list, const NTYPE list_length,
                std::vector<T2> & v)
  {
    for (NTYPE i = 0; i < list_length; i++)
      { v.push_back(list[i]); }
  }


  /// Add element to a list.  Store overflow if list is full.
  /// @param[out] overflow_list Pairs of list indices and elements
  ///    where element could not be add to list.
  template <typename ITYPE, typename ETYPE0, typename ETYPE1,
            typename NTYPE0, typename NTYPE1>
  void add_to_list_store_overflow
  (const ITYPE ilist, const ETYPE0 element, ETYPE1 * const list,
   NTYPE0 & num_elements_in_list, const NTYPE1 max_num_elements_in_list, 
   std::vector<std::pair<ITYPE,ETYPE0>> & overflow_list)
  {
    if (num_elements_in_list < max_num_elements_in_list) {
      list[num_elements_in_list] = element;
      num_elements_in_list++;
    }
    else {
      overflow_list.push_back(std::pair<ITYPE,ETYPE0>(ilist,element));
    }
  }


  // **************************************************
  // Class LIST_OF_LISTS
  // **************************************************

  /// Class representing a list of lists.
  /// Used for representing a polytope mesh, vertex-poly incidences,
  ///   vertex-edge incidences, etc.
  template <typename ETYPE, typename NTYPE>
  class LIST_OF_LISTS {

  public:

    typedef NTYPE NUM_TYPE;
    typedef ETYPE ELEMENT_TYPE;

    /// list_length[i] = Length of list i.
    std::vector<NTYPE> list_length;

    /// first_element[i] = Index in element[] of first element of list i.
    std::vector<NTYPE> first_element;

    /// List of list elements.
    std::vector<ETYPE> element;

    // constructor
    LIST_OF_LISTS(){}

    /// Number of lists.
    NTYPE NumLists() const
    { return(list_length.size()); }

    /// Number of elements of list i.
    NTYPE ListLength(const NTYPE i) const
    { return(list_length[i]); }

    /// Return index of j'th element of list ilist in vector element.
    inline NTYPE ElementIndex(const NTYPE ilist, const NTYPE j) const
    { return(first_element[ilist]+j); }

    /// Return j'th element of list ilist.
    ETYPE Element(const NTYPE ilist, const NTYPE j) const
    { return(element[ElementIndex(ilist, j)]); }

    /// Return index in element[] of first element of list ilist.
    template <typename NTYPE2>
    const NTYPE FirstElement(const NTYPE2 ilist) const
    { return(first_element[ilist]); }

    /// Return pointer to first element of list ilist.
    const ETYPE * List(const NTYPE ilist) const
    { return(&(element[first_element[ilist]])); }

    /// Return true if list ilist contains element el.
    template <typename ETYPE2>
    bool DoesListContain(const NTYPE ilist, const ETYPE2 el) const
    {
      return(does_list_contain(List(ilist), ListLength(ilist), el));
    }

    /// Return true if list ilist contains element el.
    /// - Version which returns first location of element el in list.
    /// @param[out] iloc Location of element el in list.
    template <typename ETYPE2, typename ITYPE>
    bool DoesListContain(const NTYPE ilist, const ETYPE2 el, ITYPE & iloc) const
    {
      return(does_list_contain(List(ilist), ListLength(ilist), el, iloc));
    }

    /// Return true if ilistA equals ilistB.
    /// - Lists are equal only if lists have same length 
    ///   and ilistA[k] = ilistB[k] for every k.
    bool AreListsEqual(const NTYPE ilistA, const NTYPE ilistB) const;

    /// Count number of distinct elements in list ilist.
    /// @pre List ilist is sorted.
    NTYPE CountNumDistinct(const NTYPE ilist) const;

    /// Count number of distinct elements in each list and store 
    ///   in array num_distinct[].
    /// @pre Every list is sorted.
    template <typename NTYPE2>
    void CountNumDistinct(std::vector<NTYPE2> & num_distinct) const;

    /// Copy.
    template <typename LLTYPE>
    void Copy(const LLTYPE & list_of_lists);

    /// Return non-const reference to j'th element of list ilist.
    ETYPE & ElementRef(const NTYPE ilist, const NTYPE j)
    { return(element[ElementIndex(ilist, j)]); }

    /// Return reference to const j'th element of list ilist.
    const ETYPE & ElementRefConst(const NTYPE ilist, const NTYPE j) const
    { return(element[ElementIndex(ilist, j)]); }

    /// Set number of lists.
    template <typename NTYPE2>
    void SetNumLists(const NTYPE2 num_lists);

    /// Set all list lengths to L.
    /// @pre SetNumLists() must be called before SetListLengths().
    template <typename LTYPE>
    void SetListLengths(const LTYPE L);

    /// Set length of list i to L[i].
    /// @param L[] Array of list lengths. L[i] is length of list i.
    template <typename LTYPE, typename NTYPE2>
    void SetListLengths(const LTYPE L[], const NTYPE2 num_lists);

    /// Set array first_element[].
    /// @pre Array list_length[] is set.
    void SetFirstElement();

    /// Create list of lists of uniform length.
    template <typename LTYPE, typename NTYPE2>
    void CreateUniformLengthLists(const LTYPE L, const NTYPE2 num_lists);

    /// Create list of lists.
    /// @param L[] Array of list lengths. L[i] is length of list i.
    template <typename LTYPE, typename NTYPE2>
    void CreateLists(const LTYPE L[], const NTYPE2 num_lists);

    /// Create list of lists.
    /// - Version using C++ STL vector for array L[] of list lengths.
    template <typename LTYPE>
    void CreateLists(const std::vector<LTYPE> & L);

    /// Allocate array element[].
    /// @pre Array first_element[] is set.
    void AllocArrayElement();

    /// Add a list of length list_length2.
    template <typename ETYPE2, typename NTYPE2>
    void AddList(const ETYPE2 list2_element[],
                 const NTYPE2 list2_length);

    /// Add lists where each list has length list_length2.
    /// @param num_elements Number of elements in each list.
    template <typename ETYPE2, typename NTYPE1, typename NTYPE2>
    void AddLists(const ETYPE2 list2_element[],
                  const NTYPE1 num_elements,
                  const NTYPE2 list2_length);

    /// Add lists where each list has length list_length2.
    template <typename ETYPE2>
    void AddLists(const std::vector<ETYPE2> & list2_element,
                  const NTYPE list2_length);

    /// Replace elements of each list.
    /// @param replacement[] Element i is replaced by replacement[i].
    /// @pre Each element has integer value 
    ///    in range [0..(remplacement.size()-1].    
    template <typename ETYPE2>
    void ReplaceElements(const ETYPE2 * replacement);

    /// Replace elements of each list.
    /// - C++ STL vector format for replacement.
    template <typename ETYPE2>
    void ReplaceElements(const std::vector<ETYPE2> & replacement);

    /// Sort the elements of each list.
    void SortEachList();

    /// Return list of list indices in sorted order.
    /// @param sorted_list[i] = Index of i'th list in sorted order.
    /// @pre SortEachList() should be called before GetSortedListIndices.
    template <typename ITYPE>
    void GetSortedListIndices
    (std::vector<ITYPE> & sorted_list) const;

    /// Remove all lists.
    void Clear();
  };


  // **************************************************
  // LIST compare
  // **************************************************

  /// function class for comparing lists in LIST_OF_LISTS
  /// @pre: Each list is sorted.
  template <typename LLTYPE> class LIST_LESS_THAN {
  protected:
    const LLTYPE * list_of_lists;

  public:
    LIST_LESS_THAN(const LLTYPE * list_of_lists2):
      list_of_lists(list_of_lists2)
    {};

    typedef typename LLTYPE::NUM_TYPE NUM_TYPE;
    typedef typename LLTYPE::ELEMENT_TYPE ELEMENT_TYPE;

    bool operator ()(const int i0, const int i1) const
    { 
      if (list_of_lists->list_length[i0] < list_of_lists->list_length[i1])
        { return(true); }
      else if (list_of_lists->list_length[i0] > list_of_lists->list_length[i1])
        { return(false); }
      else {
        // list_of_lists->list_length[i0] = list_of_lists->list_length[i1]
        NUM_TYPE list_length = list_of_lists->list_length[i0];

        const ELEMENT_TYPE * list0 = list_of_lists->List(i0);
        const ELEMENT_TYPE * list1 = list_of_lists->List(i1);

        return(std::lexicographical_compare
               (list0, list0+list_length, 
                list1, list1+list_length));
      }
    };
  };


  // **************************************************
  // Member functions for class LIST_OF_LISTS
  // **************************************************

  template <typename ETYPE, typename NTYPE>
  bool LIST_OF_LISTS<ETYPE,NTYPE>::
  AreListsEqual(const NTYPE ilistA, const NTYPE ilistB) const
  {
    if (ListLength(ilistA) != ListLength(ilistB)) { return(false); }
    for (NTYPE j = 0; j < ListLength(ilistA); j++) {
      if (Element(ilistA,j) != Element(ilistB,j)) 
        { return(false); }
    }

    return(true);
  }

  template <typename ETYPE, typename NTYPE>
  NTYPE LIST_OF_LISTS<ETYPE,NTYPE>::
  CountNumDistinct(const NTYPE ilist) const
  {
    NTYPE n = 0;
    if (ListLength(ilist) > 0) {
      n = 1;
      const NTYPE k = first_element[ilist];
      for (NTYPE j = 1; j < ListLength(ilist); j++) {
        if (element[k+j] != element[k+j-1])
          { n++; }
      }
    }

    return(n);
  }

  template <typename ETYPE, typename NTYPE>
  template <typename NTYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  CountNumDistinct(std::vector<NTYPE2> & num_distinct) const
  {
    
    if (num_distinct.size() != NumLists()) 
      { num_distinct.resize(NumLists()); }

    for (NTYPE i = 0; i < NumLists(); i++) 
      { num_distinct[i] = CountNumDistinct(i); }
  }


  template <typename ETYPE, typename NTYPE>
  template <typename LLTYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::Copy(const LLTYPE & list_of_lists2)
  {
    const NTYPE num_lists = list_of_lists2.NumLists();

    list_length.resize(num_lists);
    first_element.resize(num_lists);
    element.resize(list_of_lists2.element.size());

    std::copy(list_of_lists2.list_length.begin(), 
              list_of_lists2.list_length.end(),
              list_length.begin());
    std::copy(list_of_lists2.first_element.begin(), 
              list_of_lists2.first_element.end(),
              first_element.begin());
    std::copy(list_of_lists2.element.begin(), 
              list_of_lists2.element.end(),
              element.begin());
  }


  template <typename ETYPE, typename NTYPE>
  template <typename NTYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  SetNumLists(const NTYPE2 num_lists)
  {
    list_length.assign(num_lists, 0);
    first_element.resize(num_lists);
    element.clear();
  }


  template <typename ETYPE, typename NTYPE>
  template <typename LTYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  SetListLengths(const LTYPE L)
  {
    for (NTYPE i = 0; i < NumLists(); i++) 
      { list_length[i] = L; }
  }


  template <typename ETYPE, typename NTYPE>
  template <typename LTYPE, typename NTYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  SetListLengths(const LTYPE L[], const NTYPE2 num_lists)
  {
    SetNumLists(num_lists);
    for (NTYPE2 i = 0; i < NumLists(); i++) 
      { list_length[i] = L[i]; }
  }


  template <typename ETYPE, typename NTYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::SetFirstElement()
  {
    const NTYPE num_lists = list_length.size();
    if (first_element.size() != num_lists)
      { first_element.resize(num_lists); }

    if (first_element.size() == 0) { return; }

    first_element[0] = 0;
    for (NTYPE i = 1; i < num_lists; i++) {
      first_element[i] = first_element[i-1]+list_length[i-1];
    }

    element.clear();
  }

  template <typename ETYPE, typename NTYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::AllocArrayElement()
  {
    const NTYPE num_lists = first_element.size();

    if (num_lists == 0) { 
      element.clear();
      return;
    }

    element.resize(first_element.back()+list_length.back());
  }


  template <typename ETYPE, typename NTYPE>
  template <typename LTYPE, typename NTYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  CreateUniformLengthLists(const LTYPE L, const NTYPE2 num_lists)
  {
    this->SetNumLists(num_lists);
    this->SetListLengths(L);
    this->SetFirstElement();
    this->AllocArrayElement();
  }


  template <typename ETYPE, typename NTYPE>
  template <typename LTYPE, typename NTYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  CreateLists(const LTYPE L[], const NTYPE2 num_lists)
  {
    this->SetListLengths(L, num_lists);
    this->SetFirstElement();
    this->AllocArrayElement();
  }


  template <typename ETYPE, typename NTYPE>
  template <typename LTYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  CreateLists(const std::vector<LTYPE> & L)
  {
    CreateLists(IJK::vector2pointer(L), L.size());
  }


  template <typename ETYPE, typename NTYPE>
  template <typename ETYPE2, typename NTYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  AddList(const ETYPE2 list2_element[],
          const NTYPE2 list2_length)
  {
    const NTYPE ifirst = element.size();
    list_length.push_back(list2_length);
    first_element.push_back(ifirst);

    for (NTYPE j = 0; j < list2_length; j++) {
      const ETYPE el = list2_element[j];
      element.push_back(el);
    }
  }


  template <typename ETYPE, typename NTYPE>
  template <typename ETYPE2, typename NTYPE1, typename NTYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  AddLists(const ETYPE2 list2_element[],
           const NTYPE1 num_elements,
           const NTYPE2 list2_length)
  {
    const NTYPE num_lists = num_elements/list2_length;
    for (NTYPE ilist = 0; ilist < num_lists; ilist++) {
      const NTYPE ifirst = element.size();
      list_length.push_back(list2_length);
      first_element.push_back(ifirst);

      for (NTYPE j = 0; j < list2_length; j++) {
        const ETYPE el = list2_element[ilist*list2_length+j];
        element.push_back(el);
      }
    }
  }


  template <typename ETYPE, typename NTYPE>
  template <typename ETYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  AddLists(const std::vector<ETYPE2> & list2_element,
           const NTYPE list2_length)
  {
    AddLists(IJK::vector2pointer(list2_element), 
             list2_element.size(), list2_length);
  }

  template <typename ETYPE, typename NTYPE>
  template <typename ETYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  ReplaceElements(const ETYPE2 * replacement)
  {
    for (NTYPE ilist = 0; ilist < NumLists(); ilist++) {
      for (NTYPE j = 0; j < ListLength(ilist); j++) {
        const NTYPE k = ElementIndex(ilist, j);
        const ETYPE2 x = element[k];
        element[k] = replacement[x];
      }
    }
  }


  template <typename ETYPE, typename NTYPE>
  template <typename ETYPE2>
  void LIST_OF_LISTS<ETYPE,NTYPE>::
  ReplaceElements(const std::vector<ETYPE2> & replacement)
  {
    ReplaceElements(IJK::vector2pointer(replacement));
  }


  template <typename ETYPE, typename NTYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::SortEachList()
  {
    for (NTYPE ilist = 0; ilist < NumLists(); ilist++) {
      NTYPE j = first_element[ilist];

      std::sort(element.begin()+j, 
                element.begin()+j+list_length[ilist]);
    }
  }


  template <typename ETYPE, typename NTYPE>
  template <typename ITYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::GetSortedListIndices
  (std::vector<ITYPE> & sorted_list) const
  {
    LIST_LESS_THAN<LIST_OF_LISTS> list_lt(this);

    sorted_list.resize(NumLists());

    for (int i = 0; i < NumLists(); i++)
      { sorted_list[i] = i; }

    std::sort(sorted_list.begin(), sorted_list.end(), list_lt);

  }

  template <typename ETYPE, typename NTYPE>
  void LIST_OF_LISTS<ETYPE,NTYPE>::Clear()
  {
    list_length.clear();
    first_element.clear();
    element.clear();
  }


  // **************************************************
  // Functions on LIST_OF_LISTS
  // **************************************************

  /// Count number of elements with distinct tuple 0 values in list ilist.
  /// @pre ETYPE is a tuple.
  /// @pre List ilist is sorted by tuple 0 values.
  template <typename ETYPE, typename NTYPE0, typename NTYPE1>
  NTYPE1 count_num_distinct_tuple0
  (const LIST_OF_LISTS<ETYPE,NTYPE0> & list, const NTYPE1 ilist)
  {
    typedef typename std::tuple_element<0,ETYPE>::type TYPE0;

    NTYPE1 n = 0;
    if (list.ListLength(ilist) > 0) {
      n = 1;
      const NTYPE1 k = list.first_element[ilist];
      TYPE0 x0 = std::get<0>(list.element[k]);
      for (NTYPE1 j = 1; j < list.ListLength(ilist); j++) {
        TYPE0 x1 = std::get<0>(list.element[k+j]);
        if (x0 != x1) { n++; }
        x0 = x1;
      }
    }

    return(n);
  }


  /// Count number of elements with distinct tuple 0 values
  ///   in each list and store in array num_distinct[].
  /// @pre ETYPE is a tuple.
  /// @pre Each list is sorted by tuple 0 values.
  template <typename ETYPE, typename NTYPE0, typename NTYPE1>
  void count_num_distinct_tuple0
  (const LIST_OF_LISTS<ETYPE,NTYPE0> & list, 
   std::vector<NTYPE1> & num_distinct)
  {
    if (num_distinct.size() != list.NumLists()) 
      { num_distinct.resize(list.NumLists()); }

    for (NTYPE1 i = 0; i < list.NumLists(); i++) 
      { num_distinct[i] = count_num_distinct_tuple0(list, i); }
  }

}

#endif
