#ifndef _H_INTupleWriterSvc_
#define _H_INTupleWriterSvc_

#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/NTuple.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_INTupleWriterSvc("INTupleWriterSvc",  2 , 1 ); 

/*! @class INTupleWriterSvc
 @brief Proper Gaudi abstract interface class for the ntupleWriterSvc 
*/
class INTupleWriterSvc : virtual public IInterface
{  

public:

    /// setup before event processing - required for Gaudi services
    virtual StatusCode initialize ()=0;
    
    /// cleanup after event processing - required for all Gaudi services
    virtual StatusCode finalize ()=0;

    /// special version that adds a <EM>pointer</EM> to an item
    virtual StatusCode addItem(const std::string & tupleName, const std::string& itemName, const double* val)=0;

    virtual StatusCode saveNTuples()=0;

    virtual void storeRowFlag(bool flag)=0;

    virtual bool storeRowFlag()=0;

    /** store row flag by tuple Name option, retrive currrent
    @param tupleName Name of the tuple (TTree for RootTupleSvc implemetation)
    @param flag new value
    @return previous value
    If service does not implement, it is ignored (return false)
    */
    virtual bool storeRowFlag(const std::string& tupleName, bool flag)=0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_INTupleWriterSvc; }

};


#endif // _H_INTupleWriterSvc
