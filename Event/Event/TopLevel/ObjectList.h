// $Header$
#ifndef LHCBEVENT_OBJECTLIST_H
#define LHCBEVENT_OBJECTLIST_H 1


// Include files
#include <list>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ObjectContainerBase.h"
#include "tEvent/Utilities/ProcessingVersion.h"
#include "Event/TopLevel/Definitions.h"


// Externals 
extern const CLID& CLID_ObjectList;


//------------------------------------------------------------------------------
//
// ClassName:   ObjectList
//
//! This is one of the LHCb container classes. (The other is ObjectVector)

/*!           It is based on Standard Library (STL) std::list
              (see <A HREF="http://www.sgi.com/Technology/STL/">STL Programmer's Guide</A>)


              ObjectList has all functions of the std::list interface,
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
              - m_list (the STL list itself)

              Security :

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
class ObjectList : public ObjectContainerBase                                  {

public:
  typedef TYPE                                                 contained_type;
  typedef typename std::list<TYPE*>::value_type                value_type;

  typedef typename std::list<TYPE*>::reference                 reference;
  typedef typename std::list<TYPE*>::const_reference           const_reference;

  typedef typename std::list<TYPE*>::size_type                 size_type;
 
  typedef typename std::list<TYPE*>::iterator                  iterator;
  typedef typename std::list<TYPE*>::const_iterator            const_iterator;

  typedef typename std::list<TYPE*>::reverse_iterator          reverse_iterator;
  typedef typename std::list<TYPE*>::const_reverse_iterator    const_reverse_iterator;

#if defined(WIN32) || defined(_WIN32)
  typedef typename std::list<TYPE*>::_Nodeptr          link_type;
  link_type __Node(const ObjectList<TYPE>::iterator& i)       const        { return i._Mynode(); }
  link_type __Node(const ObjectList<TYPE>::const_iterator& i) const        { return i._Mynode(); }
#else
  typedef typename std::list<TYPE*>::link_type         link_type;
  link_type __Node(const ObjectList<TYPE>::iterator& i)       const             { return i.node; }
  link_type __Node(const ObjectList<TYPE>::const_iterator& i) const             { return i.node; }
#endif

public:
  /// Constructors
  ObjectList( const char* name = "ObjectList<TYPE>" )
    : ObjectContainerBase(name),
      m_list(0),
      m_processingVersion(0)                                                 { }
  ObjectList( const char* name,
              const ProcessingVersion& processing, 
              const DetectorDataObject* detectorData )
    : ObjectContainerBase( "ObjectList<TYPE>", detectorData ),
      m_list(0),
      m_processingVersion(0)                                                 { }

  /// Destructor
  virtual ~ObjectList()                                                        {
    for( ObjectList<TYPE>::iterator iter = begin(); iter != end(); iter++ ) {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container, and deleting the object
    (*iter)->setParent(0);
    delete *iter;
    }
  }

  /// Retrieve pointer to class defininition structure
  virtual const CLID& clID() const                                             { 
    return ObjectList<TYPE>::classID(); 
  }
  static const CLID& classID()                                                 {
    static CLID clid = TYPE::classID() + CLID_ObjectList; 
    return clid;
  }

  /// Clone operator
  const ObjectList<TYPE>& operator = (const ObjectList<TYPE> &right)           {
    processingVersion = right.m_processingVersion;
    detectorDataObject = right.m_detectorDataObject;
    m_list = right.m_list;
    return *this;
  }

  /// Return an iterator pointing to the beginning of the container
  ObjectList<TYPE>::iterator begin ()                                          {
    return m_list.begin();    
  }

  /// Return a const_iterator pointing to the beginning of the container
  ObjectList<TYPE>::const_iterator begin () const                              {
    return m_list.begin();
  }

  /// Return an iterator pointing to the end of the container
  ObjectList<TYPE>::iterator end ()                                            {
    return m_list.end();
  }

  /// Return a const_iterator pointing to the end of the container
  ObjectList<TYPE>::const_iterator end () const                                {
    return m_list.end();
  }

  /// Return a reverse_iterator pointing to the beginning
  ///   of the reversed container
  ObjectList<TYPE>::reverse_iterator rbegin ()                                 {
    return m_list.rbegin();
  }

  /// Return a const_reverse_iterator pointing to the beginning
  ///   of the reversed container
  ObjectList<TYPE>::const_reverse_iterator rbegin () const                     {
    return m_list.rbegin();
  }

  /// Return a reverse_iterator pointing to the end
  ///   of the reversed container
  ObjectList<TYPE>::reverse_iterator rend ()                                   {
    return m_list.rend();
  }

  /// Return a const_reverse_iterator pointing to the end
  ///   of the reversed container
  ObjectList<TYPE>::const_reverse_iterator rend () const                       {
    return m_list.rend();
  }

  /// Return the size of the container
  ///   Size means the number of objects stored in the container,
  ///     independently on the amount of information stored in each object
  ObjectList<TYPE>::size_type size () const                                    {
    return m_list.size();
  }
  /// The same as size(), return number of objects in the container
  virtual long numberOfObjects() const                                         {
    return m_list.size();
  }

  /// Return the largest possible size of the container
  ObjectList<TYPE>::size_type max_size () const                                {
    return m_list.max_size();
  }

  /// Return true if the size of the container is 0
  bool empty () const                                                          {
    return m_list.empty();
  }

  /// Return reference to the first element
  ObjectList<TYPE>::reference front ()                                         {
    return m_list.front();
  }

  /// Return const_reference to the first element
  ObjectList<TYPE>::const_reference front () const                             {
    return m_list.front();
  }

  /// Return reference to the last element
  ObjectList<TYPE>::reference back ()                                          {
    return m_list.back();
  }

  /// Return const_reference to the last element
  ObjectList<TYPE>::const_reference back () const                              {
    return m_list.back();
  }

  /// push_back = append = insert a new element at the end of the container
  void push_back( ObjectList<TYPE>::const_reference value )                    {
    if( 0 != value->parent() ) {
      const_cast<ObjectContainerBase*>(value->parent())->release(value);
    }
    value->setParent(this);
    m_list.push_back(value);
  }

  /// Add an object to the container
  virtual void add(ContainedObject* pObject)                                   {
    try {
      ObjectList<TYPE>::value_type ptr =
            dynamic_cast<ObjectList<TYPE>::value_type>(pObject);
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
    ObjectList<TYPE>::value_type position = m_list.back();
    // Set the back pointer to 0 to avoid repetitional searching
    // for the object in the container, and deleting the object
    position->setParent (0);
    delete position;
    // Removing from the container itself
    m_list.pop_back();
  }

  /// Release object from the container (the poiter will be removed
  ///   from the container, but the object itself will remain alive)
  ///     (see the method pop_back)
  virtual StatusCode release(ContainedObject* value)                           {
  // Find the object of value value
    ObjectList<TYPE>::iterator   iter;
    for( iter = begin(); iter != end(); iter++ )  {
      if( value == *iter ) {
        break;
      }
    }
    if( end() == iter )  {
      // Object cannot be released from the conatiner,
      // as it is not contained in it
      return StatusCode::FAILURE;
    }
    else  {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container and deleting the object
      (*iter)->setParent (0);
      erase(iter);
      return StatusCode::SUCCESS;
    }
  }

  /// Insert "value" before "position"
  ObjectList<TYPE>::iterator insert( ObjectList<TYPE>::iterator position,
                                     ObjectList<TYPE>::const_reference value ) {
    value->setParent(this);
    ObjectList<TYPE>::iterator i = m_list.insert(position, value);
    return i;
  }

  /// Erase the object at "position" from the container
  ///   The removed object will be deleted
  void erase( ObjectList<TYPE>::iterator position )                            {
    if( 0 != (*position)->parent() ) {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container, and deleting the object
      (*position)->setParent (0);
      delete *position;
    }
    // Removing from the container itself
    m_list.erase(position);
  }

  /// Erase the range [first, last) from the container
  ///   The removed object will be deleted
  void erase( ObjectList<TYPE>::iterator first,
              ObjectList<TYPE>::iterator last )                                {
    for( ObjectList<TYPE>::iterator iter = first; iter != last; iter++ )  {
      // Set the back pointer to 0 to avoid repetitional searching
      // for the object in the container, and deleting the object
      (*iter)->setParent (0);
      delete *iter;
    }
    // Removing from the container itself
    m_list.erase(first, last);
  }

  /// Remove objects from the container
  virtual void erase( ObjectContainerBase::iterator& start,
                      ObjectContainerBase::iterator& stop )                    {
    ObjectList<TYPE>::iterator loc_start((link_type)start.m_context);
    ObjectList<TYPE>::iterator loc_stop ((link_type)stop.m_context);
    ObjectList<TYPE>::erase(loc_start,loc_stop);
  }
/*
  /// Remove objects from "start" to "stop" from the container
  virtual void erase( ObjectContainerBase::iterator& start,
                      ObjectContainerBase::iterator& stop )                    {
    ObjectList<TYPE>::iterator loc_start = (ObjectList<TYPE>::iterator)start.m_context;
    ObjectList<TYPE>::iterator loc_stop  = (ObjectList<TYPE>::iterator)stop.m_context;
    ObjectList<TYPE>::erase(loc_start,loc_stop);
  }
*/
  /// Return distance of a given object from the beginning of its container
  ///   It correcponds to the "index" ( from 0 to size()-1 )
  ///   If "obj" not fount, return -1
  virtual long distance( const ContainedObject* obj ) const                    {
    long i = 0;
    ObjectList<TYPE>::const_iterator   iter;
    for( iter = begin(); iter != end(); iter++ )  {
      if( *iter == obj ) {
        return i;
      }
      i++;
    }
    return -1;
  }

  /// Return pointer to an object of a given distance
  virtual ContainedObject* containedObject( long dist )                        {
    long i = 0;
    ObjectList<TYPE>::const_iterator   iter;
    for( iter = begin(); iter != end(); iter++ )  {
      if( dist == i ) {
        return *iter;
      }
      i++;
    }
    return 0;
  }
  /// Return const pointer to an object of a given distance
  virtual const ContainedObject* containedObject( long dist ) const            {
    long i = 0;
    ObjectList<TYPE>::const_iterator   iter;
    for( iter = begin(); iter != end(); iter++ )  {
      if( dist == i ) {
        return *iter;
      }
      i++;
    }
    return 0;
  }

  /// Virtual functions (forwards to the concrete container definitions)

  /// Return iterator of the first element 
  virtual ObjectContainerBase::iterator first()                                {
    ObjectList<TYPE>::iterator i = m_list.begin();    
    ObjectContainerBase::iterator iter( this, *i, __Node(i) );
    return iter;
  }
  /// Return const_iterator of the first element 
  virtual ObjectContainerBase::const_iterator first() const                    {
    ObjectList<TYPE>::const_iterator i = m_list.begin();
    ObjectContainerBase::const_iterator iter( this, *i, __Node(i) );
    return iter;
  }

  /// Return iterator of the last element
  virtual ObjectContainerBase::iterator last()                                 {
    ObjectContainerBase::iterator iter( this, 0, __Node(m_list.end()) );
    return iter;
  }
  /// Return const_iterator of the last element
  virtual ObjectContainerBase::const_iterator last() const                     {
    ObjectContainerBase::const_iterator iter( this, 0, __Node(m_list.end()) );
    return iter;
  }

  /// Return iterator to the next element
  virtual ObjectContainerBase::iterator&
              next( ObjectContainerBase::iterator& prev )                      {
    ObjectList<TYPE>::iterator i((link_type)prev.m_context);
    i++;
    prev.m_context = __Node(i);
    prev.m_Ptr     = (i == m_list.end()) ? 0 : (*i);
    return prev;
  }
  /// Return const_iterator to the next element
  virtual ObjectContainerBase::const_iterator&
              next( ObjectContainerBase::const_iterator& prev ) const          {
    ObjectList<TYPE>::const_iterator i((link_type)prev.m_context);
    i++;
    prev.m_context = __Node(i);
    prev.m_Ptr     = (i == m_list.end()) ? 0 : (*i);
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
    s << "class ObjectList :    size = "
      << GlastEventField( GlastEvent::field12 )
      << size() << "\n";
    // Output the base class
    ObjectContainerBase::fillStream(s);
    s << "\n        Processing version = " << m_processingVersion;
    if ( 0 != size() ) {
      s << "\nContents of the STL list :";
      long   count = 0;
      ObjectList<TYPE>::const_iterator iter;
      for( iter = m_list.begin(); iter != m_list.end(); iter++, count++ ) {
        s << "\nIndex "
          << GlastEventField( GlastEvent::field12 )
          << count
          << " of object of type "<< **iter;
      }
    }
    return s;
  }

private:

  /// The STL list
  std::list<TYPE*>     m_list;
  /// Processing version of the container   // Cannot be in ObjectContainerBase,
  ProcessingVersion    m_processingVersion; // as it would create unfavourable
                                            // dependency of GAUDI on GlastEvent
};


#endif    // LHCBEVENT_OBJECTLIST_H
