// $Header$
#ifndef CalVec_h
#define CalVec_h

// LOCAL INLUDES
#include "CalDefs.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <vector>

namespace CalUtil {
  
  /** \brief STL vector wrapper restricts indexing to specified type.
      
  \author Zachary Fewtrell
      
  intended for use w/ index classes from CalDefs.h.

  \note default initial size of vector is 0
  \note default initial vector values are all set to 0.

  \param idx_type vector index type.
  \param val_type vector value type.
  */


  template <typename idx_type, typename val_type >
  class CalVec : private vector<val_type > {
  protected:
    typedef vector<val_type> parent_type;
    typedef size_t size_type;
    
  public:
    typedef typename parent_type::reference reference;
    typedef typename parent_type::const_reference const_reference;
    typedef typename parent_type::iterator iterator;
    typedef typename parent_type::const_iterator const_iterator;

  public:
    /// default to initial size based on idx_type, zero initialize like STL
    explicit CalVec() : parent_type(idx_type::N_VALS) {}
        
    /// just like STL, set initial size
    explicit CalVec(const size_type sz) : parent_type(sz) {}
    
    /// just like STL, set initial size & value
    CalVec(const size_type sz, const val_type &val) : parent_type(sz,val) {}
        
    /// just like STL, copy from another collection
    template <typename _iterator>
    CalVec(_iterator first, _iterator last) :
      parent_type(first,last) {}

    reference operator[] (const idx_type &idx) {
      return parent_type::operator[](idx.val());
    }
    const_reference operator[] (const idx_type &idx) const {
      return parent_type::operator[](idx.val());
    }

    reference at(const idx_type &idx) {
      return parent_type::at(idx.val());
    }
    const_reference at(const idx_type &idx) const {
      return parent_type::at(idx.val());
    }

    void resize(const size_type sz) {
      parent_type::resize(sz);
    }

    void resize(const size_type sz, const val_type &val) {
      parent_type::resize(sz, val);
    }
    
    void clear() {
      parent_type::clear();
    }
    
    size_type size() const {
      return parent_type::size();
    }
 
    const_iterator begin() const {return parent_type::begin();}
    iterator begin() {return parent_type::begin();}

    const_iterator end() const {return parent_type::end();}
    iterator end() {return parent_type::end();}
  };
};

#endif
