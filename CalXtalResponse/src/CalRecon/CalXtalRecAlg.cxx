//   $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "CalXtalRecAlg.h"
#include "XtalRecTool.h"
#include "../Xtalk/INeighborXtalkTool.h"
#include "CalXtalResponse/IXtalRecTool.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/stl_util.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Digi/CalDigi.h"


// EXTLIB INCLUDES
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD INCLUDES
#include <algorithm>

using namespace Event;
using namespace idents;
using namespace std;
using namespace CalUtil;

//static const AlgFactory<CalXtalRecAlg>  Factory;
//const IAlgFactory& CalXtalRecAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(CalXtalRecAlg);

CalXtalRecAlg::CalXtalRecAlg(const string& name, ISvcLocator* pSvcLocator):
  Algorithm(name, pSvcLocator),
  m_calDigiCol(0),
  m_calXtalRecCol(0),
  m_xtalRecTool(0),
  m_xtalkTool(0)
{
  // declare jobOptions.txt properties
  declareProperty("XtalRecToolName", m_recToolName="XtalRecTool");
  declareProperty("NeighborXtalkToolName", m_xtalkToolName="");
}

/// initialize / retrieve all needed Gaudi objects
StatusCode CalXtalRecAlg::initialize()
{
  StatusCode sc;
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  //-- JOB OPTIONS --//
  sc = setProperties();
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  //-- Xtal Recon Tool --//
  sc = toolSvc()->retrieveTool("XtalRecTool",
                               m_recToolName, 
                               m_xtalRecTool,
                               this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_recToolName << endreq;
    return sc;
  }

  //-- Neighbor Xtalk Tool --//
  if (!m_xtalkToolName.value().empty()) {
    sc = toolSvc()->retrieveTool("NeighborXtalkTool",
                                 m_xtalkToolName, 
                                 m_xtalkTool, 
                                 0); // shared by other code
    if (sc.isFailure() ) {
      msglog << MSG::ERROR << "  Unable to create " << m_xtalkToolName << endreq;
      return sc;
    }
  }

  return StatusCode::SUCCESS;
}

/**
   This function is called to do reconstruction for all cal xtals.
   for all CalDigis in one event.
   It calls retrieve() method to get access to input and output TDS data.
   It iterates over all elements in CalDigiCol and calls
   XtalRecTool::calculate() method for each xtal reconstruction.
*/
StatusCode CalXtalRecAlg::execute()
{
  StatusCode sc;

  //get access to TDS data collections
  sc = retrieve(); 
  if (sc.isFailure()) return sc;

  // non-fatal error:
  /// if there's no CalDigiCol then CalXtalRecAlg is not happening, go on w/ other algs
  if (!m_calDigiCol) return StatusCode::SUCCESS;

  // initialize neighborXtalkTool
  if (m_xtalkTool)
    m_xtalkTool->buildSignalMap(*m_calDigiCol);  
       
  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator digiIter = m_calDigiCol->begin(); 
       digiIter != m_calDigiCol->end(); digiIter++) {
    
    // if there is no digi data, then move on w/ out creating
    // recon TDS data for this xtal
    if ((*digiIter)->getReadoutCol().size() < 1) continue;
    
    const CalXtalId xtalId((*digiIter)->getPackedId());
    
    // create new object to store crystal reconstructed data  
    // use auto_ptr so it is autmatically deleted when we exit early
    // on error.
    auto_ptr<CalXtalRecData> recData(new CalXtalRecData((*digiIter)->getMode(), xtalId));
                                     
    // calculate energy in the crystal
    // used for current range only
    CalVec<FaceNum, bool> belowNoise;
    CalVec<FaceNum, bool> saturated;
    CalVec<FaceNum, bool> adcSaturated;

    // convert adc values into energy/pos
    sc = m_xtalRecTool->calculate(**digiIter,
                                  *recData,
                                  belowNoise,
                                  saturated,
                                  adcSaturated,
                                  m_xtalkTool);
    
    // single xtal may not be able to recon, is not failure condition.
    if (sc.isFailure()) continue;

    // skip any xtal w/ at least one LEX8 range below noise thold
    if (contains(belowNoise,true)) continue;
    
    // skip if recData contains no entries
    if (recData->getNReadouts() == 0) continue;

    // add new reconstructed data to the collection
    // release it from the auto_ptr so it is not deleted
    m_calXtalRecCol->push_back(recData.release());
  }

  // This is where the code to deal with the ambiguity has to go.S
  sc = dynamic_cast<XtalRecTool*>(m_xtalRecTool)->ambiguity(m_calXtalRecCol);

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
  StatusCode sc;

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
  if (sc.isFailure()) {
    delete m_calXtalRecCol;
	return sc;
  }
  
  return StatusCode::SUCCESS;
}


