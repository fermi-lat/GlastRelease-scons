// $Header$



// Include files
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"

#include "astro/GPS.h"


// GlastEvent for creating the McEvent stuff
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"

// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include <list>
#include <string>
#include <vector>
#include "GaudiKernel/ParticleProperty.h"

/*! \class TestAlg
\brief 

*/

class TestAlg : public Algorithm {

public:
    //! Constructor of this form must be provided
    TestAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();


private:
    IFlux* m_flux;
    IFluxSvc* m_fsvc; /// pointer to the flux Service 
    std::string m_source_name;
    IParticlePropertySvc * m_partSvc;
    DoubleProperty m_latitude;
    DoubleProperty m_longitude;
    StringArrayProperty m_rootplot;
};


static const AlgFactory<TestAlg>  Factory;
const IAlgFactory& CRTestAlgFactory = Factory;

//------------------------------------------------------------------------------
//
TestAlg::TestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){

    declareProperty("source_name", m_source_name="default");

}

//------------------------------------------------------------------------------
/*! */
StatusCode TestAlg::initialize() {


    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get the service
    StatusCode sc = service("FluxSvc", m_fsvc);
    m_fsvc->GPSinstance()->notifyObservers();

    return sc;
}


//------------------------------------------------------------------------------
StatusCode TestAlg::execute() {

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    

    return sc;
}


//------------------------------------------------------------------------------
StatusCode TestAlg::finalize() {

    return StatusCode::SUCCESS;
}






