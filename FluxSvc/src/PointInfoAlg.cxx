/** @file PointInfoAlg.cxx
@brief declaration and definition of the class PointInfoAlg

$Header$
Note: reverted to 1.5 by THB on 11/10/2008

*/

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/SmartDataPtr.h"

// Event for access to time
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/Exposure.h"

// to write a Tree with pointing info
#include "ntupleWriterSvc/INTupleWriterSvc.h"


//flux
#include "FluxSvc/PointingInfo.h"



/** 
* \class PointInfoAlg
*
* \brief This is an Algorithm designed to store pointing information in the tuple
* \author Toby Burnett
* 
* $Header$
*/


class PointInfoAlg : public Algorithm {
public:
    PointInfoAlg(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();


private: 
    PointingInfo m_pointing_info;

    StringProperty m_root_tree;
    BooleanProperty m_save_tuple; // set true to save

    INTupleWriterSvc* m_rootTupleSvc;

};
//------------------------------------------------------------------------


static const AlgFactory<PointInfoAlg>  Factory;
const IAlgFactory& PointInfoAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
PointInfoAlg::PointInfoAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator) 
{
    // declare properties with setProperties calls

    declareProperty("pointing_info_tree_name",  m_root_tree="MeritTuple");
    declareProperty("save_pointing_info",  m_save_tuple=false);

}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode PointInfoAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    // Use the Job options service to set the Algorithm's parameters
    setProperties();


    // get a pointer to RootTupleSvc
    if( (service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " RootTupleSvc is not available" << endreq;
        m_rootTupleSvc=0;
        sc = StatusCode::FAILURE;
    }else if( !m_root_tree.value().empty() ) {
        
        m_pointing_info.setPtTuple(m_rootTupleSvc, m_root_tree.value());
    }


    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode PointInfoAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    //
    // Purpose: set tuple items
 
    m_pointing_info.set();

   
    // put pointing stuff into the root tree
    if( m_rootTupleSvc!=0 && !m_root_tree.value().empty()){
        m_rootTupleSvc->storeRowFlag(this->m_root_tree.value(), m_save_tuple);
    }


    // Here the TDS receives the exposure data
    Event::ExposureCol* exposureDBase = new Event::ExposureCol;
    sc=eventSvc()->registerObject(EventModel::MC::ExposureCol , exposureDBase);
    if(sc.isFailure()) {
        log << MSG::ERROR << EventModel::MC::ExposureCol  
            <<" could not be entered into existing data store" << endreq;
        return sc;
    }
    exposureDBase->push_back(m_pointing_info.forTDS());

    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode PointInfoAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    return sc;
}


