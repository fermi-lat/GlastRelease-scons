#ifndef _BASECNV_H
#define _BASECNV_H 1

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"

class ISvcLocator;
class DataObject;
class IGlastDetSvc;

//extern const long int SICB_StorageType;

template <class TYPE> class CnvFactory;


/** @class BaseCnv
* @brief Base class for all GLAST Converters.
*
* Based on SICb's SicbBaseCnv
*
* $Header$
*/
class BaseCnv : public Converter      {
protected:
    // access to the GlastDetSvc
    IGlastDetSvc* m_detSvc;
    
    
public:
    
    BaseCnv(const CLID& clid, ISvcLocator* svc);
    
    ~BaseCnv() { };
    
    /// Storage type and class ID
    static long storageType()      { 
        return SICB_StorageType; 
    }
    
    virtual StatusCode initialize();
    
    virtual StatusCode finalize();
    
    /// Store TDS path to link a particular converter to an object on the TDS
    //void declareObject(const std::string& fullPath, const CLID& clid, const std::string& bank, long par=0);
    
    /// Convert the transient object to the requested representation.
    virtual StatusCode createRep(DataObject* pObject, IOpaqueAddress*& refpAddress);
    
    /// Resolve the references of the converted object. 
    virtual StatusCode fillRepRefs(IOpaqueAddress* pAddress,DataObject* pObject);
    
    /// Update the converted representation of a transient object.
    virtual StatusCode updateRep(IOpaqueAddress* pAddress, DataObject* pObject);
    
    /// Update the references of an already converted object.
    virtual StatusCode updateRepRefs(IOpaqueAddress* pAddress, DataObject* pObject);

    virtual long repSvcType() const { return SICB_StorageType; }
};


#endif // _BASECNV_H
