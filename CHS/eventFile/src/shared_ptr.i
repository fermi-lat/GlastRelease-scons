// -*- mode: c++ -*-
// This is a simplified description of the Boost shared_ptr<>
// template. It includes just the features we're interested in using
// from Python, in particular the copy constructor (which we give here
// although its left implicit in the original).
//
// We add the isNull() member function to take the place of coercion
// to boolean. There is also a SHARED_PTR() macro defined which does
// the template instantiation for you.
//
// $Header$

namespace boost {
  template<class T> class shared_ptr
  {
    
  public:
    
    shared_ptr();
    
    explicit shared_ptr(T *p);
    
    shared_ptr(const shared_ptr &p);
    
    T * operator-> () const;
    
    bool unique() const;
    
    long use_count() const;
    
    void reset();

    void swap(shared_ptr &p);

  };  // shared_ptr
}; // namespace boost

// Look, Ma! I can extend a template!
%extend boost::shared_ptr {
  bool isNull() {
    return !self->get();
  }

};

%define SHARED_PTR(pointee, pointerAlias)
// The typedef is essential; otherwise the Python object passed as
// argument will be an opaque SWIG pointer to the ahred_ptr<> proxy
// object instead of the proxy object itself.
typedef boost::shared_ptr<pointee> pointerAlias;
%template(pointerAlias) boost::shared_ptr<pointee>;
%enddef
