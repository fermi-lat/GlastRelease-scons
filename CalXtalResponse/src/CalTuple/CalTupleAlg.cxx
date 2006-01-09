// $Header$
// LOCAL INCLUDES

// GLAST INCLUDES
#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "Event/Digi/CalDigi.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"


// EXTLIB INCLUDES
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD INCLUDES
#include <cstring>

using namespace std;
using namespace CalUtil;

/** @class CalTupleAlg
    @brief populates CalTuple entry w/ info derived from CalDigis & calibrations

    CalTuple contains 1 row per GLAST event w/ the following 2 values for each xtal face
    - m_calXtalAdcPed     - pedestal subtracted adc values
    - m_calXtalFaceSignal - energy deposit estimated from single face adc value if deposit centroid is assumed to be at center of xtal.

    @author Zach Fewtrell
*/
class CalTupleAlg : public Algorithm {
public:
  CalTupleAlg(const string& name, ISvcLocator* pSvcLocator);
  /// initialize internal data members.
  StatusCode initialize();
  ///  Reconstruct ene & pos for all xtal hits in event
  StatusCode execute();
  /// required by Gaudi Algorithm class
  StatusCode finalize() {return StatusCode::SUCCESS;};
  
private:
  /// Single entry in CalTuple
  struct CalTupleEntry {
    void Clear() {memset(this,0,sizeof(CalTupleEntry));}

    int m_runId;
    int m_eventId;

    /// ped subtracted adcs
    float m_calXtalAdcPed[16][8][12][2];
        
    /// Cal Xtal Face signal in scintillated MeV units.
    float m_calXtalFaceSignal[16][8][12][2];
  };

  /// reusable store for CalTuple entries.
  CalTupleEntry m_tupleEntry;

  /// name of svc from which i get calibrations
  StringProperty m_calCalibSvcName;
  /// name of tuple to which i write
  StringProperty m_tupleName;
  /// name of file to which i write tuple
  StringProperty m_tupleFilename;

  /// pointer to CalCalibSvc object.
  ICalCalibSvc *m_calCalibSvc;  

  /// pointer to tupleWriterSvc
  INTupleWriterSvc *m_tupleWriterSvc;
};

static const AlgFactory<CalTupleAlg>  Factory;
const IAlgFactory& CalTupleAlgFactory = Factory;

CalTupleAlg::CalTupleAlg(const string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator),
  m_calCalibSvc(0),
  m_tupleWriterSvc(0)
{

  // declare jobOptions.txt properties
  declareProperty("CalCalibSvc",   m_calCalibSvcName="CalCalibSvc");
  declareProperty("tupleName",     m_tupleName="CalTuple");
  declareProperty("tupleFilename", m_tupleFilename="");
}

StatusCode CalTupleAlg::initialize() {
  StatusCode sc;
  MsgStream msglog(msgSvc(), name());
        
  //-- JOB OPTIONS --//
  sc = setProperties();
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvc." << endreq;
    return sc;
  }

  sc = service("RootTupleSvc", m_tupleWriterSvc);
  // if we can't retrieve the tuple svc just pretend we never wanted
  // it and continue anyway
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Could not locate the ntupleSvc" << endreq;
    m_tupleWriterSvc = 0;
    return StatusCode::FAILURE;
  }

  else { // tuple svc was successfully found
    // keep track of any branch creation errors
    bool branchFailure = false;
    sc = m_tupleWriterSvc->addItem(m_tupleName.value(), 
                                   "RunID", 
                                   &m_tupleEntry.m_runId,
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;
    sc = m_tupleWriterSvc->addItem(m_tupleName.value(), 
                                   "EventID", 
                                   &m_tupleEntry.m_eventId,
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;
    sc = m_tupleWriterSvc->addItem(m_tupleName.value(), 
                                   "CalXtalAdcPed[16][8][12][2]",
                                   (float*)m_tupleEntry.m_calXtalAdcPed,
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;
    sc = m_tupleWriterSvc->addItem(m_tupleName.value(), 
                                   "CalXtalFaceSignal[16][8][12][2]",
                                   (float*)m_tupleEntry.m_calXtalFaceSignal,
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;
      
    if (branchFailure) {
      msglog << MSG::ERROR << "Failure creating tuple branches" << endl;
      return StatusCode::FAILURE;
    }
  }
  
  return StatusCode::SUCCESS;
}

StatusCode CalTupleAlg::execute() {
  StatusCode sc;

  // _allways_ store a row, even if it's not empty
  m_tupleWriterSvc->storeRowFlag(true);  

  // clear tuple entry for this event
  m_tupleEntry.Clear();

  // Retrieve the Event data for this event
  SmartDataPtr<Event::EventHeader> evtHdr(eventSvc(), EventModel::EventHeader);
  if (!evtHdr) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "Failed to retrieve Event" << endreq;
  } else {
    // populate event id info
    m_tupleEntry.m_runId   = evtHdr->run();
    m_tupleEntry.m_eventId = evtHdr->event();
  }
         
  // get a pointer to the input TDS data collection
  SmartDataPtr<Event::CalDigiCol> calDigiCol(eventSvc(), EventModel::Digi::CalDigiCol);

  if (!calDigiCol) {
    if (msgSvc()->outputLevel(name()) <= MSG::VERBOSE) {
      // create msglog only when needed (for performance)
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::VERBOSE << "No CalDigi data found"
             << endreq;
    }
  }  
  else {
    // loop over all calorimeter digis in CalDigiCol
    for (Event::CalDigiCol::const_iterator digiIter = calDigiCol->begin(); 
         digiIter != calDigiCol->end(); digiIter++) {

      // get CalXtalId
      CalXtalId xtalId = (*digiIter)->getPackedId();  
      // get CalUtil::XtalIdx
      XtalIdx xtalIdx(xtalId);
      TwrNum twr = xtalIdx.getTwr();
      LyrNum lyr = xtalIdx.getLyr();
      ColNum col = xtalIdx.getCol();

      // currently allways using 1st readout
      Event::CalDigi::CalXtalReadoutCol::const_iterator ro = 
        (*digiIter)->getReadoutCol().begin();

      // PER FACE LOOP
      for (FaceNum face; face.isValid();  face++) {
        // get adc range
        RngNum rng((*ro).getRange(face)); 
        // get adc values 
        float adc = (*ro).getAdc(face);   

        // get pedestals
        // pedestals
        float ped;
        float sig, cos; // not used
        RngIdx rngIdx(xtalIdx, face, rng);
        sc = m_calCalibSvc->getPed(rngIdx, ped, sig, cos);
        if (sc.isFailure()) return sc;

        // ped subtracted ADC
        // get reference to 'real' location in big array
        float &adcPed = m_tupleEntry.m_calXtalAdcPed[twr][lyr][col][face];
        adcPed =  adc - ped; 

        // face signal
        // get reference to 'real' location in big array
        float &faceSignal  = m_tupleEntry.m_calXtalFaceSignal[twr][lyr][col][face];
        sc = m_calCalibSvc->evalFaceSignal(rngIdx,adcPed, faceSignal);
        if (sc.isFailure()) return sc;
      }
    }
  }

  return StatusCode::SUCCESS;
}

