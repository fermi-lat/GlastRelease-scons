#ifndef _LdfBaseCnv_H
#define _LdfBaseCnv_H 1

#include <vector>
#include <string>
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "ILdfCnvSvc.h"

//class IGlastDetSvc;

extern unsigned const char TEST_StorageType;

/** @class LdfBaseCnv
* @brief Base class for all GLAST EBF Converters.
*
* Based on SICb's SicbBaseCnv
*
* $Header$
*/
class LdfBaseCnv : public Converter      {
protected:
    /// Pointer to event conversion service interface
    ILdfCnvSvc*  m_CnvSvc;
    /// Leaf objects giving the locations of the objects in the data store
    std::vector<ILdfCnvSvc::Leaf>   m_leaves;
    // access to the GlastDetSvc
    //IGlastDetSvc* m_detSvc;
    
    
public:
    
    LdfBaseCnv(const CLID& clid, ISvcLocator* svc);
    
    virtual ~LdfBaseCnv() { };
    
    /// Storage type and class ID
    static const unsigned char storageType()      { 
        return TEST_StorageType; 
    }
    
    virtual StatusCode initialize();
    
    virtual StatusCode finalize();
    
    /// Store TDS path to link a particular converter to an object on the TDS
    void declareObject(const std::string& fullPath, const CLID& clid, const std::string& bank, long par=0);
    
    /// Convert the transient object to the requested representation.
    virtual StatusCode createRep(DataObject* pObject, IOpaqueAddress*& refpAddress);
    
    /// Resolve the references of the converted object. 
    virtual StatusCode fillRepRefs(IOpaqueAddress* pAddress,DataObject* pObject);
    
    /// Update the converted representation of a transient object.
    virtual StatusCode updateRep(IOpaqueAddress* pAddress, DataObject* pObject);
    
    /// Update the references of an already converted object.
    virtual StatusCode updateRepRefs(IOpaqueAddress* pAddress, DataObject* pObject);
};


#endif // _LdfBaseCnv_H
