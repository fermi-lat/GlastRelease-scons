// $Header$
#ifndef LHCBEVENT_OBJECTVECTOR_H
#define LHCBEVENT_OBJECTVECTOR_H 1


// Include files
#include <vector>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ObjectContainerBase.h"
#include "Event/Utilities/ProcessingVersion.h"
#include "Event/TopLevel/Definitions.h"


// Externals 
extern const CLID& CLID_ObjectVector;


//------------------------------------------------------------------------------
//
// ClassName:   ObjectVector
//
//

//! This is a  container class for data.

/*!     It is based on Standard Library (STL) std::vector
     (see <A HREF="http://www.sgi.com/Technology/STL/">STL Programmer's Guide</A>)


     ObjectVector has all functions of the std::vector interface,
     and in addition LHCb specific functions :
     - const CLID& classID()
     - StatusCode release(ContainedObject* value)
     - long distance( const ContainedObject* obj ) const
     - const ContainedObject* containedObject( long dist ) const
     - ContainedObject* containedObject( long dist ) const
     - long numberOfObjects() const
     - void add(ContainedObject* pObject)
     - ObjectContainerBase::const_iterator first() const
     - ObjectContainerBase::iterator first()
     - ObjectContainerBase::const_iterator last() const
     - ObjectContainerBase::iterator last()
     - const ProcessingVersion& processingVersion () const
     - void setProcessingVersion (const ProcessingVersion& value)
     - bool hasNewerVersion (const ObjectContainerBase& value) const
     - bool hasTheSameVersion (const ObjectContainerBase& value) const
     - StreamBuffer& serialize( StreamBuffer& s ) const
     - StreamBuffer& serialize( StreamBuffer& s )
     - std::ostream& fillStream( std::ostream& s ) const

     It contains these data members :
     - m_processingVersion
     - m_vector (the STL vector itself)

     Security :<br>

     Each object is allowed to belong into a single container only.
     After inserting the object into the container, it takes over
     all responsibilities for the object.  E.g. erasing the object
     from its container causes removing the object's pointer from
     the container and deleting the object itself.
 */
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


template <class TYPE>
class ObjectVector : public ObjectContainerBase                                {

public:
  typedef TYPE                                                 contained_type;
  typedef typename std::vector<TYPE*>::value_type              value_type;

  typedef typename std::vector<TYPE*>::reference               reference;
  typedef typename std::vector<TYPE*>::const_reference         const_reference;

  typedef typename std::vector<TYPE*>::size_type               size_type;

  typedef typename std::vector<TYPE*>::iterator                iterator;
  typedef typename std::vector<TYPE*>::const_iterator          const_iterator;

  typedef typename std::vector<TYPE*>::reverse_iterator        reverse_iterator;
  typedef typename std::vector<TYPE*>::const_reverse_iterator  const_reverse_iterator;

public:
  /// Constructors
  ObjectVector( const char* name = "ObjectVector<TYPE>" )
    : ObjectContainerBase(name),
      m_vector(0),
      m_processingVersion(0)                                                 { }
  ObjectVector( const ObjectVector<TYPE>& value )
    : ObjectContainerBase( "ObjectVector<TYPE>", value.detectorData ),
      m_vector(value.m_vector),
      m_processingVersion(value.m_processingVersion)                         { }
  ObjectVector( const char* name,
                const ProcessingVersion& processing,
                const DetectorDataObject* detectorData )
    : ObjectContainerBase( "ObjectVector<TYPE>", detectorData ),
      m_vector(0),
      m_processingVersion(0)                                                 { }

  /// Destructor
  virtual ~ObjectVector()                                                      {
    for( ObjectVector<TYPE>::iterator i = begin(); i != end(); i++ ) {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container, and deleting the object
      (*i)->setParent (0);
      delete *i;
    }
  }

  /// Retrieve pointer to class defininition structure
  virtual const CLID& clID() const                                             {
    return ObjectVector<TYPE>::classID();
  }
  static const CLID& classID()                                                 {
    static CLID clid = TYPE::classID() + CLID_ObjectVector;
    return clid;
  }

  /// Clone operator
  const ObjectVector<TYPE>& operator = (const ObjectVector<TYPE> &right)       {
    processingVersion = right.m_processingVersion;
    detectorDataObject = right.m_detectorDataObject;
    m_vector = right.m_vector;
    return *this;
  }

  /// Return an iterator pointing to the beginning of the container
  ObjectVector<TYPE>::iterator begin ()                                        {
    return m_vector.begin();    
  }

  /// Return a const_iterator pointing to the beginning of the container
  ObjectVector<TYPE>::const_iterator begin () const                            {
    return m_vector.begin();
  }

  /// Return an iterator pointing to the end of the container
  ObjectVector<TYPE>::iterator end ()                                          {
    return m_vector.end();
  }

  /// Return a const_iterator pointing to the end of the container
  ObjectVector<TYPE>::const_iterator end () const                              {
    return m_vector.end();
  }

  /// Return a reverse_iterator pointing to the beginning
  ///   of the reversed container
  ObjectVector<TYPE>::reverse_iterator rbegin ()                               {
    return m_vector.rbegin();
  }

  /// Return a const_reverse_iterator pointing to the beginning
  ///   of the reversed container
  ObjectVector<TYPE>::const_reverse_iterator rbegin () const                   {
    return m_vector.rbegin();
  }

  /// Return a reverse_iterator pointing to the end
  ///   of the reversed container
  ObjectVector<TYPE>::reverse_iterator rend ()                                 {
    return m_vector.rend();
  }

  /// Return a const_reverse_iterator pointing to the end
  ///   of the reversed container
  ObjectVector<TYPE>::const_reverse_iterator rend () const                     {
    return m_vector.rend();
  }

  /// Return the size of the container
  ///   Size means the number of objects stored in the container,
  ///     independently on the amount of information stored in each object
  ObjectVector<TYPE>::size_type size () const                                  {
    return m_vector.size();
  }
  /// The same as size(), return number of objects in the container
  virtual long numberOfObjects() const                                         {
    return m_vector.size();
  }

  /// Return the largest possible size of the container
  ObjectVector<TYPE>::size_type max_size () const                              {
    return m_vector.max_size();
  }

  /// Return number of elements for which memory has been allocated
  ///   It is always greater than or equal to size()
  ObjectVector<TYPE>::size_type capacity () const                              {
    return m_vector.capacity();
  }

  /// Reserve place for "value" objects in the container
  ///   If "value" is less than or equal to capacity(), this call has no effect,
  ///     otherwise, it is a request for allocation of additional memory.
  ///   If the request is successful, then capacity() is >= n,
  ///     otherwise, capacity() is unchanged.
  ///   In either case, size() is unchanged
  void reserve( ObjectVector<TYPE>::size_type value )                          {
    m_vector.reserve( value );
  }

  /// Return true if the size of the container is 0
  bool empty () const                                                          {
    return m_vector.empty();
  }

  /// Return reference to the first element
  ObjectVector<TYPE>::reference front ()                                       {
    return m_vector.front();
  }

  /// Return const_reference to the first element
  ObjectVector<TYPE>::const_reference front () const                           {
    return m_vector.front();
  }

  /// Return reference to the last element
  ObjectVector<TYPE>::reference back ()                                        {
    return m_vector.back();
  }

  /// Return const_reference to the last element
  ObjectVector<TYPE>::const_reference back () const                            {
    return m_vector.back();
  }

  /// push_back = append = insert a new element at the end of the container
  void push_back( ObjectVector<TYPE>::const_reference value )                  {
    if( 0 != value->parent() ) {
      const_cast<ObjectContainerBase*>(value->parent())->release(value);
    }
    value->setParent(this);
    m_vector.push_back(value);
  }

  /// Add an object to the container
  virtual void add(ContainedObject* pObject)                                   {
    try {
      ObjectVector<TYPE>::value_type ptr =
            dynamic_cast<ObjectVector<TYPE>::value_type>(pObject);
      if ( 0 != ptr ) {
        push_back(ptr);
      }
    }
    catch(...) {
    }
  }

  /// pop_back = remove the last element from the container
  ///   The removed object will be deleted
  ///     (see the method release)
  void pop_back ()    {
    ObjectVector<TYPE>::value_type position = m_vector.back();
    // Set the back pointer to 0 to avoid repetitional searching
    // for the object in the container, and deleting the object
    position->setParent (0);
    delete position;
    // Removing from the container itself
    m_vector.pop_back();
  }

  /// Release object from the container (the poiter will be removed
  ///   from the container, but the object itself will remain alive)
  ///     (see the method pop_back)
  virtual StatusCode release(ContainedObject* value)                           {
    // Find the object of the value value
    ObjectVector<TYPE>::iterator i;
    for( i = begin(); i != end(); i++ ) {
      if( value == *i ) {
        break;
      }
    }
    if( end() == i )  {
      // Object cannot be released from the conatiner,
      // as it is not contained in it
      return StatusCode::FAILURE;
    }
    else  {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container and deleting the object
      (*i)->setParent (0);
      erase(i);
      return StatusCode::SUCCESS;
    }
  }

  /// Insert "value" before "position"
  ObjectVector<TYPE>::iterator insert( ObjectVector<TYPE>::iterator position,
                                       ObjectVector<TYPE>::const_reference value )                 {
    value->setParent(this);
    ObjectVector<TYPE>::iterator i = m_vector.insert(position, value);
    return i;
  }

  /// Erase the object at "position" from the container
  ///   The removed object will be deleted
  void erase( ObjectVector<TYPE>::iterator position )                          {
    if( 0 != (*position)->parent() ) {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container, and deleting the object
      (*position)->setParent (0);
      delete *position;
    }
    // Removing from the container itself
    m_vector.erase(position);
  }

  /// Erase the range [first, last) from the container
  ///   The removed object will be deleted
  void erase( ObjectVector<TYPE>::iterator first,
              ObjectVector<TYPE>::iterator last )                              {
    for( ObjectVector<TYPE>::iterator i = first; i != last; i++ ) {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container, and deleting the object
      (*i)->setParent(0);
      delete *i;
    }
    // Removing from the container itself
    m_vector.erase(first, last);
  }

  /// Remove objects from "start" to "stop" from the container
  virtual void erase( ObjectContainerBase::iterator& start,
                      ObjectContainerBase::iterator& stop )                    {
    ObjectVector<TYPE>::iterator loc_start = (ObjectVector<TYPE>::iterator)start.m_context;
    ObjectVector<TYPE>::iterator loc_stop  = (ObjectVector<TYPE>::iterator)stop.m_context;
    ObjectVector<TYPE>::erase(loc_start,loc_stop);
  }

  /// Return the reference to the n'th object in the container
  ObjectVector<TYPE>::reference
          operator[] (ObjectVector<TYPE>::size_type n)                         {
    return m_vector[n];
  }

  /// Return the const_reference to the n'th object in the container
  ObjectVector<TYPE>::const_reference
          operator[] (ObjectVector<TYPE>::size_type n) const                   {
    return m_vector[n];
  }

  /// Return distance of a given object from the beginning of its container
  ///   It correcponds to the "index" ( from 0 to size()-1 )
  ///   If "obj" not fount, return -1
  virtual long distance( const ContainedObject* obj ) const                    {
    long i;
    for( i = 0; i < (long)m_vector.size(); i++ ) {
      if( m_vector[i] == obj ) {
        return i;
      }
    }
    return -1;
  }

  /// Return pointer to an object of a given distance
  virtual ContainedObject* containedObject( long dist )                        {
    return m_vector[dist];
  }
  /// Return const pointer to an object of a given distance
  virtual const ContainedObject* containedObject( long dist ) const            {
    return m_vector[dist];
  }

  /// Virtual functions (forwards to the concrete container definitions)

  /// Return iterator of the first element 
  virtual ObjectContainerBase::iterator first()                                {
    if ( m_vector.begin() == m_vector.end() )   {
      return ObjectContainerBase::iterator( this, 0, m_vector.begin() );
    }
    else  {
      return ObjectContainerBase::iterator( this, m_vector.front(),
                                                  m_vector.begin() );
    }
  }
  /// Return const_iterator of the first element 
  virtual ObjectContainerBase::const_iterator first() const                    {
    if ( m_vector.begin() == m_vector.end() )   {
      return ObjectContainerBase::const_iterator( this, 0, m_vector.begin() );
    }
    else  {
      return ObjectContainerBase::const_iterator( this, m_vector.front(),
                                                        m_vector.begin() );
    }
  }

  /// Return iterator of the last element
  virtual ObjectContainerBase::iterator last()                                 {
    return ObjectContainerBase::iterator( this, 0, m_vector.end() );
  }
  /// Return const_iterator of the last element
  virtual ObjectContainerBase::const_iterator last() const                     {
    return ObjectContainerBase::const_iterator( this, 0, m_vector.end() );
  }

  /// Return iterator to the next element
  virtual ObjectContainerBase::iterator&
          next( ObjectContainerBase::iterator& prev )                          {
    ObjectVector<TYPE>::iterator i = (ObjectVector<TYPE>::iterator)prev.m_context;
    i++;
    prev.m_Ptr     = (i == m_vector.end()) ? 0 : (*i);
    prev.m_context = i;
    return prev;
  }
  /// Return const_iterator to the next element
  virtual ObjectContainerBase::const_iterator&
          next( ObjectContainerBase::const_iterator& prev ) const              {
    ObjectVector<TYPE>::const_iterator i = (ObjectVector<TYPE>::const_iterator)prev.m_context;
    i++;
    prev.m_Ptr     = (i == m_vector.end()) ? 0 : (*i);
    prev.m_context = i;
    return prev;
  }

  /// Retrieve processing version
  const ProcessingVersion& processingVersion () const                          {
    return m_processingVersion;
  }

  /// Update processing version
  void setProcessingVersion (const ProcessingVersion& value)                   {
    m_processingVersion = value;  
  }

  /// Return true, if this container has new version, than the container "value"
  bool hasNewerVersion (const ObjectContainerBase& value) const                {
    return value.m_processingVersion < m_processingVersion;
  }

  /// Return true, if both containers have the same version
  bool hasTheSameVersion (const ObjectContainerBase& value) const              {
    return value.m_processingVersion == m_processingVersion;
  }

  /// Fill the output stream (ASCII)
  virtual std::ostream& fillStream( std::ostream& s ) const                    {
    s << "class ObjectVector :    size = "
      << GlastEventField( GlastEvent::field12 )
      << size() << "\n";
    // Output the base class
    ObjectContainerBase::fillStream(s);
    s << "\n        Processing version = " << m_processingVersion;
    if ( 0 != size() ) {
      s << "\nContents of the STL vector :";
      long   count = 0;
      ObjectVector<TYPE>::const_iterator iter;
      for( iter = m_vector.begin(); iter != m_vector.end(); iter++, count++ ) {
        s << "\nIndex "
          << GlastEventField( GlastEvent::field12 )
          << count
          << " of object of type " << **iter;
      }
    }
    return s;
  }

private:
  /// The STL vector itself
  std::vector<TYPE*>     m_vector;
  /// Processing version of the container     // Cannot be in ObjectContainerBase,
  ProcessingVersion      m_processingVersion; // as it would create unfavourable
                                              // dependency of GAUDI on GlastEvent
};


#endif    // LHCBEVENT_OBJECTVECTOR_H

