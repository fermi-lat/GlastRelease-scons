#ifndef CNV_MCEVENTCNV_H 
#define CNV_MCEVENTCNV_H 1

#include "BaseCnv.h"

class IDataDirectory;
class ISvcLocator;
class MCEvent;

// Abstract factory to create the converter
template <class T> class CnvFactory;

/** 
* @class MCEventCnv
*
* $Header$
*/
class MCEventCnv : public BaseCnv { 
    friend class CnvFactory<MCEventCnv>;
public:
    
    MCEventCnv(ISvcLocator* svc);
    
    virtual ~MCEventCnv() { };
    static const CLID& classID();
    
    /// Create object
    virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);
    
private:
};
#endif // CNV_MCEVENTCNV_H
