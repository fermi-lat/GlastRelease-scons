
#ifndef _H_ntupleWriterSvc_
#define _H_ntupleWriterSvc_


// includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"
#include "ntupleWriterSvc/INtupleWriterSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"

#include <map>

//forward declarations
template <class TYPE> class SvcFactory;

/*!  
*/
class ntupleWriterSvc :  public Service, virtual public IIncidentListener,
                        virtual public INTupleWriterSvc
{  


public:

    //------------------------------------------------------------------
    //  stuff required by a Service
    
    /// perform initializations for this service. 
    virtual StatusCode initialize ();
    
    /// 
    virtual StatusCode finalize ();
    //------------------------------------------------------------------
    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    
 
   /// Query interface
   virtual StatusCode queryInterface( const IID& riid, void** ppvUnknown );

   /// Provide the named ntuple Ptr from the data store
   virtual SmartDataPtr<NTuple::Tuple> getNTuple(const char *tupleName);

   /// add a new item to an ntuple
   virtual StatusCode addItem(const char *tupleName, const char *item, const float val, bool first);

protected: 

    /// Standard Constructor
    ntupleWriterSvc ( const std::string& name, ISvcLocator* al );
    

private:

    INTupleSvc *ntupleSvc;

    // Allow SvcFactory to instantiate the service.
    friend class SvcFactory<ntupleWriterSvc>;
    /// routine that is called when we reach the end of an event
    void endEvent();
    /// create the ntuple
    StatusCode bookNTuple(int index);
    /// write an event's data into the ntuple
    StatusCode writeNTuple(int index);
    /// add an entry into an ntuple, where the item has already been declared 
    StatusCode addValue(const char *tupleName, const char *item, const float val);


    ///Maintain a list of each item so that this algorithm could handle multiple ntuples

    /// set non-blank to enable tuples
    std::vector<std::string> m_tuple_name;
    /// of the form "/1"
    std::vector<std::string> m_TDS_tuple_name;  
    /// path into the data store where the ntuple is located
    std::vector<std::string> m_fileName;

    /// Store ntuple Id and ntuple path in data store
    std::map<std::string, std::string> m_tuples;

    static unsigned int m_tupleCounter;
};


#endif // _H_ntupleWriterSvc