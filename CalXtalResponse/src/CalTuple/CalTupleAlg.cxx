// $Header$

/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "../Xtalk/INeighborXtalkTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST INCLUDES
#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "Event/Digi/CalDigi.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD INCLUDES
#include <cstring>


using namespace std;
using namespace CalUtil;
using namespace Event;
using namespace CalibData;
using namespace idents;

/** @class CalTupleAlg
    @brief populates CalTuple entry w/ info derived from CalDigis & calibrations

    CalTuple contains 1 row per GLAST event w/ the following 2 values for each xtal face
    - m_calXtalAdcRng     - adc range for matching AdcPed val
    - m_calXtalFaceSignal - energy deposit estimated from single face adc value if deposit centroid is assumed to be at center of xtal. (BESTRANGE)
    - m_calXtalFaceSignalAllRange - face signal for all 4 ranges.
    - m_calXtalAdcPed     - pedestal subtracted adc values.  (BEASTRANGE)
    - m_calXtalAdcPedAllRange - ped subtraced adc values (all available ranges, indexed by rangeId)

    jobOptions:
    - CalCalibSvc           - Cal calibration source (default="CalCalibSvc")
    - tupleName             - TTree object name (default="CalTuple")
    - tupleFilename         - Optional tuple filename (default is combined with other RootTupleSvc output)
    - NeighborXtalkToolName - (optional) used to correct for bias caused by neighboring crystal electronic crosstalk. (default="" - disabled)

    @author Z.Fewtrell
*/
class CalTupleAlg : public Algorithm {
public:
  CalTupleAlg(const string& name, ISvcLocator* pSvcLocator);
  /// initialize internal data members.
  StatusCode initialize();
  ///  Reconstruct ene & pos for all xtal hits in event
  StatusCode execute();
  /// required by Gaudi Algorithm class
  StatusCode finalize() {return StatusCode::SUCCESS;}
  
private:
  /// Single entry in CalTuple
  struct CalTupleEntry {
    CalTupleEntry() {
      Clear();
    }

    ~CalTupleEntry() {}

    void Clear() {
      m_runId = 0;
      m_eventId = 0;

      memset(m_calXtalAdcRng,     0, sizeof(m_calXtalAdcRng));
      memset(m_calXtalFaceSignal, 0, sizeof(m_calXtalFaceSignal));
      memset(m_calXtalFaceSignalAllRange, 0, sizeof(m_calXtalFaceSignalAllRange));
      memset(m_calXtalAdcPed,     0, sizeof(m_calXtalAdcPed));
      memset(m_calXtalAdcPedAllRange, 0, sizeof(m_calXtalAdcPedAllRange));

    }

    int m_runId;
    int m_eventId;

    /// \brief ped subtracted adcs
    /// 
    /// stored as a pointer b/c array size changes based on 
    float m_calXtalAdcPed[16][8][12][2];
    float m_calXtalAdcPedAllRange[16][8][12][2][4];

    /// adc range selection
    int m_calXtalAdcRng[16][8][12][2];
        
    /// Cal Xtal Face signal in scintillated MeV units. (BESTRANGE)
    float m_calXtalFaceSignal[16][8][12][2];
    /// Cal Xtal Face signal in scintillated MeV units. (ALLRANGE)
    float m_calXtalFaceSignalAllRange[16][8][12][2][4];
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

  /// pointer to NeighborXtalk tool for calculating electronic crosstalk from neighboring crystals
  INeighborXtalkTool *m_xtalkTool;

  /// name of INeigbhorXtalkTool instantiation
  StringProperty m_xtalkToolName;
};

//static const AlgFactory<CalTupleAlg>  Factory;
//const IAlgFactory& CalTupleAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(CalTupleAlg);

CalTupleAlg::CalTupleAlg(const string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator),
  m_calCalibSvc(0),
  m_tupleWriterSvc(0),
  m_xtalkTool(0)
{

  // declare jobOptions.txt properties
  declareProperty("CalCalibSvc",   m_calCalibSvcName="CalCalibSvc");
  declareProperty("tupleName",     m_tupleName="CalTuple");
  declareProperty("tupleFilename", m_tupleFilename="");
  declareProperty("NeighborXtalkToolName", m_xtalkToolName="");
}

/// init / retreive all Gaudi objects
StatusCode CalTupleAlg::initialize() {
  StatusCode sc;
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;
        
  //-- JOB OPTIONS --//
  sc = setProperties();
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName << endreq;
    return sc;
  }

  sc = service("RootTupleSvc", m_tupleWriterSvc);
  // if we can't retrieve the tuple svc just pretend we never wanted
  // it and continue anyway
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Could not locate the ntupleSvc" << endreq;
    m_tupleWriterSvc = 0;
    return sc;
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
                                   "CalXtalAdcPedAllRange[16][8][12][2][4]",
                                   (float*)m_tupleEntry.m_calXtalAdcPedAllRange,
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;

    sc = m_tupleWriterSvc->addItem(m_tupleName.value(), 
                                   "CalXtalAdcRng[16][8][12][2]",
                                   (int*)(&m_tupleEntry.m_calXtalAdcRng),
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;


    sc = m_tupleWriterSvc->addItem(m_tupleName.value(), 
                                   "CalXtalFaceSignal[16][8][12][2]",
                                   (float*)m_tupleEntry.m_calXtalFaceSignal,
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;

    sc = m_tupleWriterSvc->addItem(m_tupleName.value(), 
                                   "CalXtalFaceSignalAllRange[16][8][12][2][4]",
                                   (float*)m_tupleEntry.m_calXtalFaceSignalAllRange,
                                   m_tupleFilename);
    if (sc.isFailure()) branchFailure |= true;
      
    if (branchFailure) {
      msglog << MSG::ERROR << "Failure creating tuple branches" << endl;
      return StatusCode::FAILURE;
    }
  }

  //-- Neighbor Xtalk Tool --//
  if (!m_xtalkToolName.value().empty()) {
    sc = toolSvc()->retrieveTool(m_xtalkToolName, 
                                 m_xtalkTool,
                                 0); // shared by other code
    if (sc.isFailure() ) {
      msglog << MSG::ERROR << "  Unable to create " << m_xtalkToolName << endreq;
      return sc;
    }
  }
  
  return StatusCode::SUCCESS;
}

/** Loop through each digi
    - store adc for each digi in tuple
    - store adc range in tuple
    - calculate faceSiganl (from CalCalibSvc) 
    - optionally apply crosstalk
    - store face signal in tuple
 */
StatusCode CalTupleAlg::execute() {
  StatusCode sc;

  // _allways_ store a row, even if it's not empty
  m_tupleWriterSvc->storeRowFlag(m_tupleName.value(),
                                 true);  

  // clear tuple entry for this event
  m_tupleEntry.Clear();


  // Retrieve the Event data for this event
  SmartDataPtr<EventHeader> evtHdr(eventSvc(), EventModel::EventHeader);
  if (!evtHdr) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "Failed to retrieve Event" << endreq;
  } else {
    // populate event id info
    m_tupleEntry.m_runId   = evtHdr->run();
    m_tupleEntry.m_eventId = evtHdr->event();
  }
         
  // get a pointer to the input TDS data collection
  SmartDataPtr<CalDigiCol> calDigiCol(eventSvc(), EventModel::Digi::CalDigiCol);

  if (!calDigiCol) {
    if (msgSvc()->outputLevel(name()) <= MSG::VERBOSE) {
      // create msglog only when needed (for performance)
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::VERBOSE << "No CalDigi data found"
             << endreq;
    }
  }  
  else {
    // initialize neighborXtalkTool
    if (m_xtalkTool)
      m_xtalkTool->buildSignalMap(*calDigiCol);  

          
    // loop over all calorimeter digis in CalDigiCol
    for (CalDigiCol::const_iterator digiIter = calDigiCol->begin(); 
         digiIter != calDigiCol->end(); digiIter++) {

      // get CalXtalId
      const CalXtalId xtalId = (*digiIter)->getPackedId();  
      // get CalUtil::XtalIdx
      const XtalIdx xtalIdx(xtalId);
      const TwrNum twr = xtalIdx.getTwr();
      const LyrNum lyr = xtalIdx.getLyr();
      const ColNum col = xtalIdx.getCol();

      // usually using only 1st readout
      CalDigi::CalXtalReadoutCol::const_iterator ro = 
        (*digiIter)->getReadoutCol().begin();

      // PER FACE LOOP
      for (FaceNum face; face.isValid();  face++) {
        // get adc range
        const RngNum rng((*ro).getRange(face)); 
        // get adc values 
        const float adc = (*ro).getAdc(face);   
        const RngIdx rngIdx(xtalIdx, face, rng);

        // adc range
        m_tupleEntry.m_calXtalAdcRng[twr.val()][lyr.val()][col.val()][face.val()] = rng.val();

        // get pedestals
        // pedestals
	float ped;
        sc = m_calCalibSvc->getPed(rngIdx,ped);
        if (sc.isFailure()) return StatusCode::FAILURE;

        // ped subtracted ADC
        const float adcPed = adc - ped;

        //-- face signal --//
        if (adcPed > 0) {
          float faceSignal;
          sc = m_calCalibSvc->evalFaceSignal(rngIdx, adcPed, faceSignal);
          if (sc.isFailure()) return sc;
                  
          // (optional) Neighbor Xtalk Correction
          if (m_xtalkTool) {
            float xtalkMeV;
            sc = m_xtalkTool->calcXtalkMeV(DiodeIdx(xtalIdx,face,rng.getDiode()), xtalkMeV);
            if (sc.isFailure()) return sc;

            faceSignal -= xtalkMeV;
          }
                  
          faceSignal = max<float>(0,faceSignal);
          m_tupleEntry.m_calXtalFaceSignal[twr.val()][lyr.val()][col.val()][face.val()] = faceSignal;
          m_tupleEntry.m_calXtalFaceSignalAllRange[twr.val()][lyr.val()][col.val()][face.val()][rng.val()] = faceSignal;
        }


        // fill in 1st readout for both bestrange and allrange arrays
        m_tupleEntry.m_calXtalAdcPed[twr.val()][lyr.val()][col.val()][face.val()] = adcPed;
        m_tupleEntry.m_calXtalAdcPedAllRange[twr.val()][lyr.val()][col.val()][face.val()][rng.val()] = adcPed;

        // loop through remaining 3 readouts
        for (unsigned char nRO = 1; nRO < 4; nRO++) {
          // create my own local vars, sos i don't mess up the 
          // ones used only for 1st readout
          const CalDigi::CalXtalReadout *ro = (*digiIter)->getXtalReadout(nRO);
          if (!ro) continue;

          const RngNum rng(ro->getRange(face));
          const short adc = ro->getAdc(face);
            
          // get pedestals
          // pedestals
          const RngIdx rngIdx(xtalIdx, face, rng);
	  float ped;
          sc = m_calCalibSvc->getPed(rngIdx,ped);
          if (sc.isFailure()) return StatusCode::FAILURE;
            
          // ped subtracted ADC
          const float adcPed = adc - ped;
            
          m_tupleEntry.m_calXtalAdcPedAllRange[twr.val()][lyr.val()][col.val()][face.val()][rng.val()] = adcPed;

          //-- FACE SIGNAL --//
          if (adcPed > 0) {
            float faceSignal; 
            sc = m_calCalibSvc->evalFaceSignal(rngIdx, adcPed, faceSignal);
            if (sc.isFailure()) return sc;

            // (optional) Neighbor Xtalk Correction
            if (m_xtalkTool) {
              float xtalkMeV;
              sc = m_xtalkTool->calcXtalkMeV(DiodeIdx(xtalIdx,face,rng.getDiode()), xtalkMeV);
              if (sc.isFailure()) return sc;

              faceSignal -= xtalkMeV;
            }

            faceSignal = max<float>(0,faceSignal);
            m_tupleEntry.m_calXtalFaceSignalAllRange[twr.val()][lyr.val()][col.val()][face.val()][rng.val()] = faceSignal;
          }
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

    
