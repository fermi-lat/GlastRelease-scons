#ifndef _H_INTupleWriterSvc_
#define _H_INTupleWriterSvc_

// includes
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/NTuple.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_INTupleWriterSvc(1400, 1 , 0); 

/*! Proper Gaudi abstract interface class for the ntupleWriterSvc 
*/
class INTupleWriterSvc : virtual public IInterface
{  

public:

    virtual StatusCode initialize ()=0;
    
    virtual StatusCode finalize ()=0;

    virtual SmartDataPtr<NTuple::Tuple> getNTuple(const char *tupleName)=0;

    virtual StatusCode addItem(const char *tupleName, const char *item, const float val)=0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_INTupleWriterSvc; }

};


#endif // _H_INTupleWriterSvc
