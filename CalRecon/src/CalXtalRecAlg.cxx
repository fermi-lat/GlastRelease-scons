// File and version Information:
//   $Header$


// LOCAL INCLUDES
#include "CalXtalRecAlg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"
#include "Event/TopLevel/EventModel.h"

// EXTLIB INCLUDES
#include "CLHEP/Geometry/Transform3D.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD INCLUDES

using namespace Event;
using namespace idents;
using namespace std;

static const AlgFactory<CalXtalRecAlg>  Factory;
const IAlgFactory& CalXtalRecAlgFactory = Factory;

CalXtalRecAlg::CalXtalRecAlg(const string& name, ISvcLocator* pSvcLocator):
  Algorithm(name, pSvcLocator)
{
  declareProperty("xtalEneToolName", m_eneToolName="XtalEneTool");
  declareProperty("xtalPosToolName", m_posToolName="XtalPosTool");
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

  //-- EXTRACT INT CONSTANTS --//
  double value;
  // map containing pointers to integer constants to be read
  // with their symbolic names from xml file used as a key 
  typedef map<int*,string> PARAMAP;
  PARAMAP param;

  //     filling the map with information on constants to be read 
  param[&m_xNum]         = string("xNum");
  param[&m_yNum]         = string("yNum");
  param[&m_eTowerCAL]    = string("eTowerCAL");
  param[&m_eLATTowers]   = string("eLATTowers");
  param[&m_nCsISeg]      = string("nCsISeg");
  param[&m_eXtal]        = string("eXtal");
    
  //-- RETRIEVE GlastDevSvc --//

  sc = service("GlastDetSvc", m_detSvc);
  // loop over all constants information contained in the map
  for(PARAMAP::iterator iter=param.begin(); iter!=param.end();iter++){
    //  retrieve constant
    if(!m_detSvc->getNumericConstByName((*iter).second, &value)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*iter).second
             <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*iter).first)= int(value); // store retrieved value 
  }
        
  //-- EXTRACT DOUBLE CONSTANTS --//

  // map containing pointers to double constants to be read
  // with their symbolic names from xml file used as a key 
  typedef map<double*,string> DPARAMAP;
  DPARAMAP dparam; 

  dparam[&m_CsILength]  = string("CsILength");
    
  for(DPARAMAP::iterator dIter=dparam.begin(); dIter!=dparam.end();dIter++){
    if(!m_detSvc->getNumericConstByName((*dIter).second,(*dIter).first)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*dIter).second << " not defined" << endreq;
      return StatusCode::FAILURE;
    } 
  }

  //-- JOB OPTIONS --//
  sc = setProperties();
  if (sc.isFailure()) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  //-- CalXtalResponse TOOLS --//
  sc = toolSvc()->retrieveTool(m_eneToolName,m_xtalEneTool);
  if (sc.isFailure() ) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "  Unable to create " << m_eneToolName << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool(m_posToolName,m_xtalPosTool);
  if (sc.isFailure() ) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "  Unable to create " << m_posToolName << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}


/**
   This function is called to do reconstruction
   for all hitted crystals in one event.
   It calls retrieve() method to get access to input and output TDS data.
   It iterates over all elements in CalDigiCol and calls
   computeEnergy() and computePosition() methods doing real reconstruction.
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
  
  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator digiIter = m_calDigiCol->begin(); 
       digiIter != m_calDigiCol->end(); digiIter++) {

    // if there is no digi data, then move on w/ out creating
    // recon TDS data for this xtal
    if ((*digiIter)->getReadoutCol().size() < 1) continue;

    CalXtalId xtalId = (*digiIter)->getPackedId();

    // create new object to store crystal reconstructed data     
    CalXtalRecData* recData = 
      new CalXtalRecData((*digiIter)->getMode(),xtalId);
           
    // calculate energy in the crystal
    bool below_thresh;
    sc = computeEnergy(*recData, **digiIter, below_thresh);
    if (sc.isFailure()) return sc;

    if(!below_thresh){      
      // calculate position in the crystal
      sc = computePosition(*recData, **digiIter);   
      if (sc.isFailure()) return sc;

      // add new reconstructed data to the collection
      m_calXtalRecCol->push_back(recData);
    } else {
      delete recData;
    }
  }

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

  // get a pointer to the input TDS data collection
  m_calDigiCol = SmartDataPtr<CalDigiCol>(eventSvc(),
                                          EventModel::Digi::CalDigiCol);
  if (!m_calDigiCol) {
    if (msgSvc()->outputLevel(name()) <= MSG::VERBOSE) {
      // create msglog only when needed for performance
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::VERBOSE << "No CalDigi data found"
             << endreq;
    }
  }

  m_calXtalRecCol = 0;

  //  create output data collection
  m_calXtalRecCol = new CalXtalRecCol;

  DataObject* pnode=0;

  // search for CalRecon section of Event directory in TDS
  sc = eventSvc()->retrieveObject( EventModel::CalRecon::Event, pnode );
    
  // if the required directory doesn't exist - create it
  if( sc.isFailure() ) {
    sc = eventSvc()->registerObject( EventModel::CalRecon::Event,
                                     new DataObject);
    if( sc.isFailure() ) {
      // if cannot create the directory - write an error message
      // create msglog only when needed for performance
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << "Could not create CalRecon directory"
             << endreq;
      return sc;
    }
  }
    
  //register output data collection as a TDS object
  sc = eventSvc()->registerObject(EventModel::CalRecon::CalXtalRecCol,
                                  m_calXtalRecCol);
  if (sc.isFailure()) return sc;
  
  return StatusCode::SUCCESS;
}

StatusCode CalXtalRecAlg::computeEnergy(CalXtalRecData &recData, 
                                        const Event::CalDigi &digi,
                                        bool &below_thresh)
{
  StatusCode sc;

  CalXtalId xtalId = digi.getPackedId();

  const CalDigi::CalXtalReadoutCol& roCol = digi.getReadoutCol();
  
  // xtal wide threshold flag
  below_thresh = false;    

  // currently allways using 1st readout
  CalDigi::CalXtalReadoutCol::const_iterator ro = roCol.begin();

  // get readout range number for both crystal faces
  CalXtalId::AdcRange rangeP = 
    (CalXtalId::AdcRange)(*ro).getRange(CalXtalId::POS); 
  CalXtalId::AdcRange rangeM = 
    (CalXtalId::AdcRange)(*ro).getRange(CalXtalId::NEG); 
  
  // get adc values 
  int adcP = (*ro).getAdc(CalXtalId::POS);   
  int adcM = (*ro).getAdc(CalXtalId::NEG);   

  float ene;

  // used for current range only
  bool range_below_thresh = false;
  if (msgSvc()->outputLevel(name()) <= MSG::VERBOSE) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::VERBOSE;
    msglog.stream() << "id=" << xtalId // needed to support setw() manipulator
                    << "\trangeP=" << int(rangeP) << " adcP=" << setw(4) <<adcP 
                    << "\trangeN=" << int(rangeM) << " adcM=" << setw(4) <<adcM;
    msglog << endreq;
  } 

  // convert adc values into energy
  sc = m_xtalEneTool->calculate(xtalId,
                                rangeP, rangeM,
                                adcP,  adcM,
                                ene, 
                                range_below_thresh,
                                below_thresh);
  if (sc.isFailure()) return sc;

  if (range_below_thresh) {
    below_thresh = true;
    return StatusCode::SUCCESS;
  }
    
  // create output object
  CalXtalRecData::CalRangeRecData* rangeRec =
    new CalXtalRecData::CalRangeRecData(rangeP,ene,rangeM,ene);

  if (msgSvc()->outputLevel(name()) <= MSG::VERBOSE) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::VERBOSE
           << " xtalId=" << xtalId
           << " ene=" << ene
           << endreq;
  }
  
  // add output object to output collection
  recData.addRangeRecData(*rangeRec);
  delete rangeRec;
  
  return StatusCode::SUCCESS;
}



StatusCode CalXtalRecAlg::computePosition(CalXtalRecData &recData, 
                                          const CalDigi &digi)
{
  MsgStream msg(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;

  // get crystal identification 
  CalXtalId xtalId = digi.getPackedId();

  //-- CONSTRUCT GEOMETRY VECTORS FOR XTAL --//
  // unpack crystal identification into tower, layer and column number
  int layer = xtalId.getLayer();
  int tower = xtalId.getTower();
  int col   = xtalId.getColumn();

  // create Volume Identifier for segment 0 of this crystal
  idents::VolumeIdentifier segm0Id;
  segm0Id.append(m_eLATTowers);
  segm0Id.append(tower/m_xNum);
  segm0Id.append(tower%m_xNum);
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(layer);
  segm0Id.append(layer%2); 
  segm0Id.append(col);
  segm0Id.append(m_eXtal);
  segm0Id.append(0);

  HepTransform3D transf;

  //get 3D transformation for segment 0 of this crystal
  m_detSvc->getTransform3DByID(segm0Id,&transf);
  //get position of the center of the segment 0
  Vector vect0 = transf.getTranslation();
  // create Volume Identifier for the last segment of this crystal
  idents::VolumeIdentifier segm11Id;
  // copy all fields from segm0Id, except segment number
  for(int ifield = 0; ifield<CalDefs::fSegment; ifield++)
    segm11Id.append(segm0Id[ifield]);
  segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
  //get 3D transformation for the last segment of this crystal
  m_detSvc->getTransform3DByID(segm11Id,&transf);
  //get position of the center of the last segment
  Vector vect11 = transf.getTranslation();

  Point p0(0.,0.,0.);           
  // position of the crystal center
  Point pCenter = p0+(vect0+vect11)*0.5; 
  //normalized vector of the crystal direction 
  Vector dirXtal = (vect11-vect0)*m_nCsISeg/(m_nCsISeg-1);  

  //-- GET BEST RANGE READOUT --//
  // readout 0 will contain the best ranges in 1-range mode
  // and ALL-RANGE mode
  // theoretically it wont work in 4-range readout mode w/ out
  // range selection, but that mode is rarely used.
  const CalDigi::CalXtalReadoutCol& roCol = digi.getReadoutCol();
  const CalDigi::CalXtalReadout &ro = *roCol.begin();

  int adcP = ro.getAdc(CalXtalId::POS);
  int adcN = ro.getAdc(CalXtalId::NEG);
  CalXtalId::AdcRange rngP = (CalXtalId::AdcRange)ro.getRange(CalXtalId::POS);
  CalXtalId::AdcRange rngN = (CalXtalId::AdcRange)ro.getRange(CalXtalId::NEG);

  float pos;
  sc = m_xtalPosTool->calculate(xtalId,rngP,rngN,adcP,adcN,pos);
  if (sc.isFailure()) return sc;

  if (msgSvc()->outputLevel(name()) <= MSG::VERBOSE) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::VERBOSE
           << " xtalId=" << xtalId
           << " pos=" << pos
           << endreq;
  }

  // put 1D position info into 3D vector
  // 'pos' is in units of xtal Length, convert to rel units (-1->1)
  pos /= m_CsILength; 
  Point pXtal = pCenter+dirXtal*pos;

  // store calculated position in the reconstructed data
  // for the best readout range
  (recData.getRangeRecData(0))->setPosition(pXtal);

  return StatusCode::SUCCESS;
}
