// $Header$
#define ITEMCNV_CPP 


// Include files
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/Interfaces/IDataDirectory.h"
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/DataSvc/RegistryEntry.h"
#include "Address.h"
#include "ItemCnv.h"


//------------------------------------------------------------------------------
//
// Implementation of class :  ItemCnv<TYPE>
//
// Author :                   Markus Frank
//
//------------------------------------------------------------------------------


/// Create the transient representation of an object.
/*
template <class TYPE>
StatusCode ItemCnv<TYPE>::createObj(IOpaqueAddress* pAddress,
                                          DataObject*& refpObject)   {
  TYPE* new_object  = new TYPE();
  //StatusCode status = updateObj(pAddress, new_object);
  //if ( StatusCode::SUCCESS == status )    {
    refpObject = new_object;
    return StatusCode::SUCCESS;
    //return status;
  //}
  //delete new_object;
  //return status;
}
*/

/// Resolve the references of the created transient object.
//template <class TYPE>
//StatusCode ItemCnv<TYPE>::fillObjRefs(IOpaqueAddress* pAddress, DataObject* pObject)    {
//  return updateObjRefs( pAddress, pObject );
//}


/// Update the transient object from the other representation.
/*
template <class TYPE>
StatusCode ItemCnv<TYPE>::updateObj(IOpaqueAddress* pAddress, 
                                          DataObject* pObject)   {
  // This is where the Monte Carlo particles in SICB were loaded in one big go...
    return StatusCode::SUCCESS;
}
*/

/// Update the references of an updated transient object.
//template <class TYPE>
//StatusCode ItemCnv<TYPE>::updateObjRefs(IOpaqueAddress* pAddress, 
//                                              DataObject* pObject)    {
//  StatusCode status   = StatusCode::SUCCESS;
//  return StatusCode::SUCCESS;  // just return success...
//}


/// Convert the transient object to the requested representation.
//template <class TYPE>
//StatusCode ItemCnv<TYPE>::createRep(DataObject* pObject, IOpaqueAddress*& refpAddress)
//{
//  return StatusCode::FAILURE;
//}


/// Resolve the references of the converted object. 
//template <class TYPE>
//StatusCode ItemCnv<TYPE>::fillRepRefs(IOpaqueAddress* pAddress, DataObject* pObject)
//{
//  return StatusCode::FAILURE;
//}


/// Update the converted representation of a transient object.
//template <class TYPE>
//StatusCode ItemCnv<TYPE>::updateRep(IOpaqueAddress* pAddress, DataObject* pObject)
//{
 // return StatusCode::FAILURE;
//}


/// Update the references of an already converted object.
//template <class TYPE>
//StatusCode ItemCnv<TYPE>::updateRepRefs(IOpaqueAddress* pAddress, DataObject* pObject)
//{
//  return StatusCode::FAILURE;//
//}


/// Helper: Add data leaf
/*
template <class TYPE>
StatusCode ItemCnv<TYPE>::addDataLeaf( RegistryEntry* base, 
                                           RegistryEntry* ent,
                                           IOpaqueAddress* addr )   {
    // HMA This routine is not called since our Event object is simple
}
*/

/// Initialize the ItemCnv<TYPE>
/*
template <class TYPE>
StatusCode ItemCnv<TYPE>::initialize()                          {
  // Put in here Converter initialization
  return Converter::initialize();
}
*/

/// Standard Constructor
/*
template <class TYPE>
ItemCnv<TYPE>::ItemCnv(ISvcLocator* svc)
: BaseCnv(TYPE::classID(), svc)   { }
*/
  //: Converter(SICB_StorageType, TYPE::classID(), svc)               { }
