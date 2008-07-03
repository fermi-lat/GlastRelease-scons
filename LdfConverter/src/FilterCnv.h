#ifndef CNV_FilterCnv_H 
#define CNV_FilterCnv_H 1

#include "LdfBaseCnv.h"

class IDataDirectory;
class ISvcLocator;
namespace Event {
}
// Abstract factory to create the converter
template <class T> class CnvFactory;

/** @class FilterCnv
* @brief Concrete converter for the Filter header on the TDS /Event/Filter
*
* $Header$
*/
class FilterCnv : public LdfBaseCnv { 
    friend class CnvFactory<FilterCnv>;
public:
    
    FilterCnv(ISvcLocator* svc);
    
    virtual ~FilterCnv() { };
    static const CLID& classID();
    
    /// Create object
    virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);
    
private:
};
#endif // CNV_FilterCnv_H

