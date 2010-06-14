// $Header$

/*
 * HISTORY
 *
 * DATE         WHO    WHAT
 * 5/28/03  russell    Added what I deemed to be a kludge to take care of the 
 *                     fact that the algorithm gets called make after all the
 *                     events have been read. If any of the DIGI pointers are
 *                     NULL, I take this as bad, and just return -1. I think
 *                     getting called back when there are no more events to
 *                     process is a bug.
 * 6/11/03  navid      Renamed this package to EbfWriter
 */
           

#include <stdio.h>

// Gaudi system includes
#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/Digi/TkrDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/AcdDigi.h"

#include "Event/Filter/Ebf.h"

#include "EbfAcdData.h"
#include "EbfCalData.h"
#include "EbfTkrData.h"
#include "EbfGltData.h"
#include "EbfGltCounters.h"
#include "EbfOutput.h"


#include <fstream>
#include <cassert>

/** 
 * @class EbfWriter
 * @brief An algorithm to convert the digi data to ebf
 * 
 * $Header$
*/
class EbfWriter : public Algorithm 
{
public:
    EbfWriter(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    /// read in the Event Header information
    void printHeader(const Event::EventHeader& );

    /// read the MCEvent header information
    void printMcHeader(const Event::MCEvent&);

    /// read the Monte Carlo Particle data
    void printMcParticles(const Event::McParticleCol&);

    /// read the Monte Carlo Position Hit data
    void printPositionHits(const Event::McPositionHitCol&);

    /// read the Monte Carlo Integrating Hit data
    void printIntegratingHits(const Event::McIntegratingHitCol&);

    /// read the hit strip data
    void printTkrDigi(const Event::TkrDigiCol&);

    /// read in CAL ADC
    void printCalDigi(const Event::CalDigiCol&);

    /// read in ACD ADC 
    void printAcdDigi(const Event::AcdDigiCol&);

    EbfOutput m_output;

    EbfGltCounters  m_latcounters;
    

    /// parameter to store the maximum size of an event
    /// this should be fairly static
    int         m_maxEvtSize;
    
    EbfCalConstants m_calConstants;
    std::ofstream output;
};



static const AlgFactory<EbfWriter>  Factory;
const IAlgFactory& EbfWriterFactory = Factory;



EbfWriter::EbfWriter(const std::string& name, ISvcLocator* pSvcLocator):Algorithm(name, pSvcLocator)
{
    // declare properties with setProperties calls
    declareProperty("MaxEvtSize", m_maxEvtSize=0x10000);

    return;
    
}



//! set parameters and attach to various perhaps useful services.
StatusCode EbfWriter::initialize()
{

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream      log(msgSvc(), name());

    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    IGlastDetSvc *detSvc;
    sc = service ("GlastDetSvc", detSvc);
    if (sc.isFailure()) return sc;
    m_calConstants.initialize (detSvc, log);

    m_output.setPrint (log.level () <= MSG::DEBUG);
    m_output.open ("", m_maxEvtSize);

    /* Set up the EBF time stuff */
    m_latcounters.initialize ((unsigned int)2*20000000, // T0 = 2 seconds
                              (unsigned int)0,   // Initialize dead time value
                              (double)20000000., // LAT clock frequency
                              (double).00001,    // Amp of variation,
                                                 // 1 part in 10**-5
                              (double)90 * 60.); // Period of the clock
                                                 // variation, in secs
   return sc;
}



StatusCode EbfWriter::execute()
{
    //
    // Procedure and Method:  
    //     Process one event, extracting all the colections from the TDS, and
    //     passing them (if they exist) to individual functions
    //


    MsgStream  log (msgSvc(), name());
    MSG::Level level = log.level ();
    int        print = level <= MSG::INFO;
    StatusCode    sc = StatusCode::SUCCESS;    
    
    DataObject *pnode=0;
    if(eventSvc()->retrieveObject(EventModel::Filter::Event,pnode).isFailure())
      eventSvc()->registerObject(EventModel::Filter::Event,new DataObject);

    Event::Ebf *newEbf=new Event::Ebf;
    eventSvc()->registerObject(EventModel::Filter::Ebf, newEbf);
    
    //
    // TKR
    // ---
    //
    SmartDataPtr<Event::TkrDigiCol> 
                 tkrDigiCol(eventSvc(), EventModel::Digi::TkrDigiCol);

    // If no tracker data... 
    if (!tkrDigiCol) return StatusCode::FAILURE;
    
    EbfTkrData tkr;
    tkr.fill  (tkrDigiCol);


    //
    // CAL
    // ---
    //
    SmartDataPtr<Event::CalDigiCol> 
                 calDigiCol(eventSvc(), EventModel::Digi::CalDigiCol);

    // If no CAL data... 
    if (!calDigiCol) return StatusCode::FAILURE;

    EbfCalData   cal;        
    cal.fill    (calDigiCol, m_calConstants);


    //
    // ACD
    // ---
    //
    SmartDataPtr<Event::AcdDigiCol> 
                 acdDigiCol(eventSvc(), EventModel::Digi::AcdDigiCol);

    // If no ACD data... 
    if (!acdDigiCol) return StatusCode::FAILURE;

    EbfAcdData acd;    
    acd.fill  (acdDigiCol);


    //
    // GLT
    // ---
    // 
    SmartDataPtr<Event::EventHeader>header(eventSvc(), 
                                           EventModel::EventHeader);
    if (!header) return StatusCode::FAILURE;

    EbfGltData glt;
    glt.fill (header->event(), 
              header->time(), 
              &m_latcounters,
              &acd, 
              &tkr,
              &cal);


    //
    // Put the contributor's data into EBF format and write it out 
    //
    m_output.format (&acd, &cal, &tkr, &glt);
    if (level <= MSG::INFO) m_output.print  ();
    unsigned int length;
    char *data=m_output.write  (length);
    newEbf->set(data,length);
    return sc;
}


/// clean up
StatusCode EbfWriter::finalize()
{
    return StatusCode::SUCCESS;
}
