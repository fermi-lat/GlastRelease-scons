// $Header$
#ifndef ITEMCNV_H
#define ITEMCNV_H 1


// Include files
#include "BaseCnv.h"


// Externals
extern unsigned const char SICB_StorageType;


//------------------------------------------------------------------------------
//
// ClassName :    ItemCnv<TYPE>
//  
// Description :  Converter of a single item
//
// Modified SicbItemCnv
//
//------------------------------------------------------------------------------


template <class TYPE>
class ItemCnv : public BaseCnv { //public Converter	    {

public:

  /// Create the transient representation of an object.
 virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject) {
    TYPE* new_object  = new TYPE();
  //StatusCode status = updateObj(pAddress, new_object);
  //if ( StatusCode::SUCCESS == status )    {
    refpObject = new_object;
    return StatusCode::SUCCESS;
  };

  /// Resolve the references of the created transient object.
  //virtual StatusCode fillObjRefs(IOpaqueAddress* pAddress, DataObject* pObject);

  /// Update the transient object from the other representation.
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* refpObject) {
      return StatusCode::SUCCESS;
  }

  /// Update the references of an updated transient object.
  //virtual StatusCode updateObjRefs(IOpaqueAddress* pAddress, DataObject* pObject);

  /// Convert the transient object to the requested representation.
  //virtual StatusCode createRep(DataObject* pObject, IOpaqueAddress*& refpAddress);

  /// Resolve the references of the converted object. 
  //virtual StatusCode fillRepRefs(IOpaqueAddress* pAddress,DataObject* pObject);

  /// Update the converted representation of a transient object.
  //virtual StatusCode updateRep(IOpaqueAddress* pAddress, DataObject* pObject);

  /// Update the references of an already converted object.
  //virtual StatusCode updateRepRefs(IOpaqueAddress* pAddress, DataObject* pObject);

  /// class ID
  static const CLID& classID()                  { return TYPE::classID(); }

protected:

  /// Standard Constructor
  ItemCnv(ISvcLocator* svc)
      : BaseCnv(TYPE::classID(), svc) {};

  /// Standard Destructor
  virtual ~ItemCnv()                        { }

  virtual StatusCode updateObj(int* data, TYPE* pObject)  {
    return StatusCode::SUCCESS;
  }

  /// Helper: Add data leaf
  //StatusCode addDataLeaf (RegistryEntry* base, RegistryEntry* ent, IOpaqueAddress* addr);

//  virtual StatusCode updateObj(int* data, TYPE* pObject) = 0;

  //virtual StatusCode updateObjRefs(int* data, TYPE* pObject) = 0;

};


#endif // ITEMCNV_H 
