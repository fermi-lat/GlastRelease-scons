// $Header$
#ifndef Address_H
#define Address_H 1


// Include files
#include <string>
#include <vector>
#include "Gaudi/Kernel/StatusCode.h"
#include "Gaudi/Kernel/GenericAddress.h"


// Forward declarations
class IDataDirectory;


// Externals
extern unsigned char SICB_StorageType;


//------------------------------------------------------------------------------
//
// ClassName :    Address
//  
// Description :  Definition of a dummy address
//
//
//------------------------------------------------------------------------------


class Address : public GenericAddress   {
public:


  /// Map of references contained in this object
  class Entry  {
  public:
    int   m_id;
    int*  m_ptr;
    void* m_obj;
    Entry(int id, int* ptr, void* obj) : m_id(id), m_ptr(ptr), m_obj(obj)   {
    }
    int id()  {
      return m_id;
    }
    int* pointer()  {
      return m_ptr;
    }
    template <class TYPE> void loadPointer(TYPE*& pObject)  {
      pObject = static_cast<TYPE*>(m_obj);
    }
  };
  typedef std::vector<Entry>    References;

protected:
  References      m_refs;
  int m_recid;  //THB: does this make sense?
public:
  /// Validate address
  StatusCode validate();
  /// Standard Constructor
  Address(const CLID& cl, const std::string& fname, int recid, const std::string& bank);
  /// Constructor for non-initialized addresses
  Address(const GenericLink& link);
  /// Constructor for non-initialized addresses
  Address(const GenericLink* link);
  /// Standard destructor
  virtual ~Address()    {
  }
};


#endif  // Address_H
