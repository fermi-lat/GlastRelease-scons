/** @file SetSeeds.cxx
@brief declaration and definition of the class SetSeeds

$Header$

*/

// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

// Event for creating the McEvent stuff
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastRandomSvc.h"

/** 
* \class SetSeeds
*
* \brief This is an Algorithm designed to set the event seeds by accessing a function member in GlastRandomSvc 
* \author Toby Burnett
* 
* $Header$
*/

class SetSeeds : public Algorithm {
public:
    SetSeeds(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    
};
//------------------------------------------------------------------------

//static const AlgFactory<SetSeeds>  Factory;
//const IAlgFactory& SetSeedsFactory = Factory;
DECLARE_ALGORITHM_FACTORY( SetSeeds );


//------------------------------------------------------------------------
//! ctor
SetSeeds::SetSeeds(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator) 
{
    // declare properties with setProperties calls
}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode SetSeeds::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
   
    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode SetSeeds::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    SmartDataPtr<Event::MCEvent> mcheader(eventSvc(), EventModel::MC::Event);
    //if (mcheader == 0) {
    if (!mcheader) {
            log << MSG::ERROR << " no MC header" << endreq;
            return StatusCode::FAILURE;
    }
    log << MSG::DEBUG;
    if( log.isActive()) mcheader->fillStream(log.stream());
    log << endreq;

    Event::McParticleCol*  pcol=  
        SmartDataPtr<Event::McParticleCol>(eventSvc(), "/Event/MC/McParticleCol");

    if( pcol==0){ 
        log<< MSG::ERROR << "No source of particles!" << endreq;
        return StatusCode::FAILURE;
    }  
    int particle_count = pcol->size();  
    if( particle_count ==1){

        log << MSG::DEBUG;
        if( log.isActive()){
            Event::McParticleCol::const_iterator piter=pcol->begin();
            (*piter)->fillStream(log.stream());
        }
        log << endreq;
    }else{
        log<< MSG::ERROR << "Expecting a single particle: found " << particle_count << endreq;
    }

    
    // get the event header to set the time
    Event::EventHeader* h = 0; 

    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    //if(0==header) {
    if(!header) {
        // not already there: try to register instead
        sc = eventSvc()->registerObject(EventModel::EventHeader, h=new Event::EventHeader);
        if( sc.isFailure()) {
            log << MSG::WARNING << " could not find or register the event header" << endreq;
        }
    } else{  h = header;
    }

        // get run and seq from mcheader
    int run=mcheader->getRunNumber(), 
        seq=mcheader->getSequence();

   // set up random seeds
    GlastRandomSvc::instance()->applySeeds( run, seq);

    // copy stuff from the MC header to the Event header
    h->setRun(run);
    h->setEvent(seq );
    h->setTime(mcheader->time());
    h->setTrigger(0); // flag that not already set

    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode SetSeeds::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    return sc;
}

