// File and version Information:
//   $Header$

// LOCAL INCLUDES
#include "CalXtalRecAlg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"
#include "Event/TopLevel/EventModel.h"

// EXTLIB INCLUDES
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD INCLUDES

using namespace Event;
using namespace idents;
using namespace std;
using namespace CalDefs;

static const AlgFactory<CalXtalRecAlg>  Factory;
const IAlgFactory& CalXtalRecAlgFactory = Factory;

const char *CalTupleEntry::m_tupleDesc = "RunID/i:EventID/i:CalXtalAdcPed[16][8][12][2]/F:CalXtalFaceSignal[16][8][12][2]/F";


CalXtalRecAlg::CalXtalRecAlg(const string& name, ISvcLocator* pSvcLocator):
  Algorithm(name, pSvcLocator),
  m_calDigiCol(0),
  m_calXtalRecCol(0),
  m_xtalRecTool(0),
  m_tupleFile(0),
  m_tupleBranch(0),
  m_tupleTree(0)
{
  declareProperty("XtalRecToolName", m_recToolName="XtalRecTool");
  declareProperty("tupleFilename",   m_tupleFilename="");
}

/** 
    This function sets values to private data members,
    representing the calorimeter geometry and digitization
    constants. Information  from xml files is obtained using 
    GlastdetSvc::getNumericConstByName() function.
    To make this constant extraction in a loop, the pairs
    'constant pointer, constant name' are stored in
    map container. <p>
    Double and integer constants are extracted separatly,
    because constants of both types are returned
    by getNumericConstByName() as double.
*/      
StatusCode CalXtalRecAlg::initialize()
{
  StatusCode sc;
  MsgStream msglog(msgSvc(), name());
        
  //-- JOB OPTIONS --//
  sc = setProperties();
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  //-- Xtal Recon Tool --//
  sc = toolSvc()->retrieveTool(m_recToolName, m_xtalRecTool, this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_recToolName << endreq;
    return sc;
  }

  //-- Open Tuple File (Optional) --//
  if (m_tupleFilename.value() != "") {
    m_tupleFile = new TFile(m_tupleFilename.value().c_str(),"RECREATE","CalTuple");
    if (!m_tupleFile) {
      msglog << MSG::ERROR << "Unable to create TFile object: " << m_tupleFilename << endreq;
      return StatusCode::FAILURE;
    }
         
    m_tupleTree = new TTree("CalTuple","CalTuple");
    if (!m_tupleTree) {
      msglog << MSG::ERROR << "Unable to create CalTuple TTree object" << endreq;
      return StatusCode::FAILURE;
    }
         
    m_tupleBranch = m_tupleTree->Branch("CalTupleEntry",      // branch name
                                        &m_tupleEntry,        // fill object
                                        m_tupleEntry.m_tupleDesc);
    if (!m_tupleBranch) { 
      // dunno how big to make it so i made it small, enough for one
      // object + any overhead.
         
      msglog << MSG::ERROR << "Couldn't create tuple branch" << endreq;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

/**
   This function is called to do reconstruction
   for all hitted crystals in one event.
   It calls retrieve() method to get access to input and output TDS data.
   It iterates over all elements in CalDigiCol and calls
   computeXtal() and computePosition() methods doing real reconstruction.
*/
StatusCode CalXtalRecAlg::execute()
{
  StatusCode sc = StatusCode::SUCCESS;

  //get access to TDS data collections
  sc = retrieve(); 
  // non-fatal error:
  /// if there's no CalDigiCol then CalXtalRecAlg is not happening, go on w/ other algs
  if (!m_calDigiCol) return StatusCode::SUCCESS;
  // fatal error  if (sc.isFailure()) return sc;

  // reset optional tuple entry
  if (m_tupleTree) {
    m_tupleEntry.Clear();
    
    m_tupleEntry.m_runId   = m_evtHdr->run();
    m_tupleEntry.m_eventId = m_evtHdr->event();
  }
       
  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator digiIter = m_calDigiCol->begin(); 
       digiIter != m_calDigiCol->end(); digiIter++) {
    
    // if there is no digi data, then move on w/ out creating
    // recon TDS data for this xtal
    if ((*digiIter)->getReadoutCol().size() < 1) continue;
    
    CalXtalId xtalId = (*digiIter)->getPackedId();
    
    // create new object to store crystal reconstructed data  
    // use auto_ptr so it is autmatically deleted when we exit early
    // on error.
    auto_ptr<CalXtalRecData> recData(new CalXtalRecData((*digiIter)->getMode(), xtalId));
                                     
    // calculate energy in the crystal
    // used for current range only
    bool belowThreshP    = false;
    bool belowThreshM    = false;
    bool xtalBelowThresh = false;
    bool saturatedP      = false;
    bool saturatedM      = false;

    // convert adc values into energy/pos
    sc = m_xtalRecTool->calculate(*m_evtHdr,
                                  **digiIter,
                                  *recData,
                                  belowThreshP,
                                  belowThreshM,
                                  xtalBelowThresh,
                                  saturatedP,
                                  saturatedM,
                                  (m_tupleFile) ? &m_tupleEntry : 0 // optional tuple entry
                                  );
    // single xtal may not be able to recon, is not failure condition.
    if (sc.isFailure()) continue;

    // skip any xtal w/ at least one LEX8 range below threshold
    if (xtalBelowThresh) continue;

    // add new reconstructed data to the collection
    // release it from the auto_ptr so it is not deleted
    m_calXtalRecCol->push_back(recData.release());
  }

  // fill optional tuple
  if (m_tupleTree)
    m_tupleTree->Fill();

  return StatusCode::SUCCESS;
}

/** 
    Purpose and method: 
    This function provides access to the TDS input and output data
    by setting the private data members m_calDigiCol            
    and m_calXtalRecCol
    
    TDS input: CalDigiCol
    TDS output: CalXtalrecCol
*/
StatusCode CalXtalRecAlg::retrieve()
{
  StatusCode sc = StatusCode::SUCCESS;


  // Retrieve the Event data for this event
  m_evtHdr = SmartDataPtr<Event::EventHeader>(eventSvc(), EventModel::EventHeader);
  if (!m_evtHdr) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "Failed to retrieve Event" << endreq;
    return StatusCode::FAILURE;
  }
         
  // get a pointer to the input TDS data collection
  m_calDigiCol = SmartDataPtr<CalDigiCol>(eventSvc(),
                                          EventModel::Digi::CalDigiCol);
  if (!m_calDigiCol) {
    if (msgSvc()->outputLevel(name()) <= MSG::VERBOSE) {
      // create msglog only when needed (for performance)
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::VERBOSE << "No CalDigi data found"
             << endreq;
    }
  }

  m_calXtalRecCol = 0;

  DataObject* pnode=0;

  // search for CalRecon section of Event directory in TDS
  sc = eventSvc()->retrieveObject( EventModel::CalRecon::Event, pnode );
    
  // if the required directory doesn't exist - create it
  if( sc.isFailure() ) {
    sc = eventSvc()->registerObject( EventModel::CalRecon::Event, new DataObject);
    if( sc.isFailure() ) {
      // if cannot create the directory - write an error message
      // create msglog only when needed (for performance)
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << "Could not create CalRecon directory"
             << endreq;
      return sc;
    }
  }
    
  //  create output data collection
  m_calXtalRecCol = new CalXtalRecCol;

  //register output data collection as a TDS object
  sc = eventSvc()->registerObject(EventModel::CalRecon::CalXtalRecCol,
                                  m_calXtalRecCol);
  if (sc.isFailure()) return sc;

  
  return StatusCode::SUCCESS;
}


StatusCode CalXtalRecAlg::finalize() {    
  // make sure optional tuple is closed out
  if (m_tupleFile) {
    m_tupleFile->Write();
    m_tupleFile->Close(); // trees deleted
    delete m_tupleFile;
  }

  if (m_xtalRecTool)
    m_xtalRecTool->finalize();

  return StatusCode::SUCCESS;
}
