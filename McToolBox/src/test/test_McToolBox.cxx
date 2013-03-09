// $Header$
// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

// TDS class declarations: input data, and McParticle tree

#include "Event/TopLevel/EventModel.h"

#include "Event/Digi/TkrDigi.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"



// Define the class here instead of in a header file: 
//  not needed anywhere but here!
//----------------------------------------------------
/** 
* test_McToolBox
*
* @brief  A miminal test of McToolBox, using as few other packages as possible
*
* @author Leon Rochester
*
* $Header$
*/

class test_McToolBox : public Algorithm {
public:
    test_McToolBox(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    //! number of times called
    int m_count; 
    //! the GlastDetSvc used for access to detector info
};
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_McToolBox );

//static const AlgFactory<test_McToolBox>  Factory;
//const IAlgFactory& test_McToolBoxFactory = Factory;
DECLARE_ALGORITHM_FACTORY(test_McToolBox);

//------------------------------------------------------------------------
//! ctor
test_McToolBox::test_McToolBox(const std::string& name, 
                             ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_McToolBox::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_McToolBox::execute()
{
    
    // First stab a a test program
    // can be fleshed out as required
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << endreq <<  "Call " << ++m_count << ": " ;
        
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_McToolBox::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    //log  << MSG::INFO << m_count << " call(s)." << endreq;
    
    return sc;
}



