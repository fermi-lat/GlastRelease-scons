// -*- mode: c++ -*-
/**
 * @file eventFile.in
 * @brief Interface file for SWIG-generated wrapper.
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Id$
 */

%module py_eventFile

%{
#include <stdexcept>
#include <iostream>
#include "boost/shared_ptr.hpp"
#include "boost/utility.hpp"
#include "eventFile/LSEReader.h"
#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/LPA_Handler.h"
#include "eventFile/EBF_Data.h"
#include "eventFile/LSE_Keys.h"
#include "LSE_Event.h"

using namespace eventFile;
%}

// use STL types
%include stl.i
%include std_vector.i

// gate exceptions back to the interpreter
%exception {
  try {
    $action
  } catch ( std::exception& eObj ) {
    PyErr_SetString( PyExc_RuntimeError, const_cast<char*>(eObj.what()) );
    return NULL;
  }
};

// make a template for some method returns
%template(UnsignedVector) std::vector< unsigned >;
%template(UnsignedPair) std::pair< unsigned, unsigned >;

// wrap the event-structure components
%include ../eventFile/LPA_Handler.h
%template(HandlerVector) std::vector< LPA_Handler >;

%include ../eventFile/LSE_GemTime.h
%include ../eventFile/LSE_Context.h
%include ../eventFile/LSE_Info.h
%include ../eventFile/EBF_Data.h
%include ../eventFile/LSE_Keys.h

// wrap the event-container object 
%include LSE_Event.h

// wrap the smart-pointer to the event container
%include shared_ptr.i
SHARED_PTR(LSE_Event, EventPtr);

// wrap the file-reader object, exposing only what we need
%newobject eventFile::LSEReader::nextEvent;
class eventFile::LSEReader : private boost::noncopyable {
public:
  LSEReader( const std::string& );
  virtual ~LSEReader();

  /// wrap all the accessors
  unsigned runid() const;
  unsigned begSec() const;
  unsigned endSec() const;
  unsigned long long evtcnt() const;
  unsigned long long begGEM() const;
  unsigned long long endGEM() const;
  std::pair<unsigned, unsigned> seqErr( int ) const;
  std::pair<unsigned, unsigned> dfiErr( int ) const;

};

// extend the file-reader
%extend eventFile::LSEReader {
  /// get a container of all the data for the next event
  const EventPtr nextEvent()
  {
    eventFile::LSE_Context        ctx;
    eventFile::EBF_Data           ebf;
    eventFile::LSE_Info::InfoType itype;
    eventFile::LPA_Info           pinfo;
    eventFile::LCI_ACD_Info       ainfo;
    eventFile::LCI_CAL_Info       cinfo;
    eventFile::LCI_TKR_Info       tinfo;
    eventFile::LSE_Keys::KeysType ktype;
    eventFile::LPA_Keys           pakeys;
    eventFile::LCI_Keys           cikeys;
    if ( self->read( ctx, ebf, itype, pinfo, ainfo, cinfo, tinfo, ktype, pakeys, cikeys ) ) {
      return EventPtr( new LSE_Event( ctx, ebf, itype, pinfo, ainfo, cinfo, tinfo, ktype, pakeys, cikeys ) );
    } else {
      return EventPtr();
    }
  }
};

// python-friendly access to the EBF event data
%extend eventFile::EBF_Data {
  // Create a python buffer object that refers to the internal
  // buffer containing the EBF data
  PyObject* getData() {
    return PyBuffer_FromMemory( const_cast< unsigned char* >( self->data() ), self->size() );
  }
  // Create a Python string containing a copy of the EBF data
  PyObject* copyData() {
    return PyString_FromStringAndSize( reinterpret_cast< const char* >( self->data() ), self->size() );
  }
};
