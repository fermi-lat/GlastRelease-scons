// File and version Information:
//   $Header$
//
// Description:
//    CalXtalRecAlg is an algorithm to reconstruct calorimeter
//    information in each individual crystal
//
// Author: A.Chekhtman

#include "CalXtalRecAlg.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/Digi/CalDigi.h"
#include "Event/TopLevel/EventModel.h"
#include "CLHEP/Geometry/Transform3D.h"

static const AlgFactory<CalXtalRecAlg>  Factory;
const IAlgFactory& CalXtalRecAlgFactory = Factory;
using namespace Event;
using namespace idents;
using namespace std;

// constructor
CalXtalRecAlg::CalXtalRecAlg(const string& name, ISvcLocator* pSvcLocator):
  Algorithm(name, pSvcLocator)
{
  declareProperty("xtalEneToolName", m_eneToolName="XtalEneTool");
  declareProperty("xtalPosToolName", m_posToolName="XtalPosTool");
}


StatusCode CalXtalRecAlg::initialize()

  // Purpose and method:
  //           This function sets values to private data members,
  //           representing the calorimeter geometry and digitization
  //           constants. Information  from xml files is obtained using 
  //           GlastdetSvc::getNumericConstByName() function.
  //           To make this constant extraction in a loop, the pairs
  //           'constant pointer, constant name' are stored in
  //           map container. 
  //           Double and integer constants are extracted separatly,
  //           because constants of both types are returned
  //           by getNumericConstByName() as double.
  //      
        
        
{
  MsgStream log(msgSvc(), name());
  StatusCode sc;

  // extracting int constants
  double value;  // intermediate variable for reading constants from
  // xml file as doubles and converting them to interger 
  typedef map<int*,string> PARAMAP;
  PARAMAP param; // map containing pointers to integer constants to be read
  // with their symbolic names from xml file used as a key 

  //     filling the map with information on constants to be read 
    
  param[&m_xNum]         = string("xNum");
  param[&m_yNum]         = string("yNum");
  param[&m_eTowerCAL]    = string("eTowerCAL");
  param[&m_eLATTowers]   = string("eLATTowers");
  param[&m_CALnLayer]    = string("CALnLayer");
  param[&m_nCsIPerLayer] = string("nCsIPerLayer");
  param[&m_nCsISeg]      = string("nCsISeg");
  param[&m_eXtal]        = string("eXtal");
    
  // now try to find the GlastDevSvc service
    
  //    IGlastDetSvc* detSvc;
  sc = service("GlastDetSvc", m_detSvc);
  // loop over all constants information contained in the map
  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){
    //  attempt to get the constant value via the method of GlastDetSvc
    if(!m_detSvc->getNumericConstByName((*it).second, &value)) {
      // if not successful - give the error message and return
      log << MSG::ERROR << " constant " <<(*it).second
          <<" not defined" << endreq;
      return StatusCode::FAILURE;
      //  if successful - fill the constant using the pointer from the map
    } else *((*it).first)= int(value);
  }
        
  // extracting double constants
  typedef map<double*,string> DPARAMAP;
  DPARAMAP dparam; // map containing pointers to double constants to be read
  // with their symbolic names from xml file used as a key 

  dparam[&m_CsILength]  = string("CsILength");
    
  for(DPARAMAP::iterator dit=dparam.begin(); dit!=dparam.end();dit++){
    if(!m_detSvc->getNumericConstByName((*dit).second,(*dit).first)) {
      log << MSG::ERROR << " constant " <<(*dit).second << " not defined" << endreq;
      return StatusCode::FAILURE;
    } 
  }

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool(m_eneToolName,m_xtalEneTool);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "  Unable to create " << m_eneToolName << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool(m_posToolName,m_xtalPosTool);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "  Unable to create " << m_posToolName << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}


StatusCode CalXtalRecAlg::execute()

  //  Purpose and method:
  //     This function is called to do reconstruction
  //     for all hitted crystals in one event.
  //     It calls retrive() method to get access to input and output TDS data.
  //     It iterates over all elements in CalDigiCol and calls
  //     computeEnergy() and computePosition() methods doing real reconstruction.
  //
  //   TDS input: 
  //            CalDigiCol* m_calDigiCol - private class member containing
  //            a pointer to the calorimeter digi callection 
  //
  //   TDS output: 
  //            CalXtalRecCol* m_calXtalRecCol - private class member containing
  //            a pointer to the calorimeter crystal reconstructed data collection 
  //
  //     

{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());

  sc = retrieve(); //get access to TDS data collections

        
  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator it = m_calDigiCol->begin(); 
       it != m_calDigiCol->end(); it++) {
    CalXtalId xtalId = (*it)->getPackedId();
                
    // create new object to store crystal reconstructed data     
    CalXtalRecData* recData = 
      new CalXtalRecData((*it)->getMode(),xtalId);
           
           
    // calculate energy in the crystal
    bool below_thresh;
    sc = computeEnergy(recData, *it, below_thresh);
    if (sc.isFailure()) return sc;

    if(!below_thresh){      
      // calculate position in the crystal
      computePosition(recData, *it);   

      // add new reconstructed data to the collection
      m_calXtalRecCol->push_back(recData);
    } else {
      delete recData;
    }
  }

  return sc;
}


StatusCode CalXtalRecAlg::finalize()
  // empty function: required by base class (Algorithm)
{
  StatusCode sc = StatusCode::SUCCESS;
  return sc;
}



StatusCode CalXtalRecAlg::retrieve()

  // Purpose and method: 
  //            This function provides access to the TDS input and output data
  //            by setting the private data members m_calDigiCol            
  //            and m_calXtalRecCol
  //
  //   TDS input: CalDigiCol
  //   TDS output: CalXtalrecCol
  //        
  //        
        
{
        
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;

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
      log << MSG::ERROR << "Could not create CalRecon directory"
          << endreq;
      return sc;
    }
  }

    
  // get a pointer to the input TDS data collection
  m_calDigiCol = SmartDataPtr<CalDigiCol>(eventSvc(),
                                          EventModel::Digi::CalDigiCol); 

  //register output data collection as a TDS object
  sc = eventSvc()->registerObject(EventModel::CalRecon::CalXtalRecCol,
                                  m_calXtalRecCol);
  return sc;
}




StatusCode CalXtalRecAlg::computeEnergy(CalXtalRecData* recData, const CalDigi* digi, bool &below_thresh)

  //   Purpose and method:
  //                 This function calculates the energy for one crystal.
  //                 It makes a loop over all readout ranges (1 or 4, depending
  //                 on readout mode), and converts adc values for both crystal
  //                 faces into energy, using constants contained in private data
  //
  //  
  //  Input: CalDigi* digi - pointer to digitized calorimeter data for
  //                         one crystal
  //                            
  //  Output: CalXtalRecData* recData - pointer to reconstructed data for
  //                                    this crystal
                                      



{
  MsgStream log(msgSvc(), name());
  StatusCode sc;

  CalXtalId xtalId = digi->getPackedId();

  const CalDigi::CalXtalReadoutCol& roCol = digi->getReadoutCol();
  
  below_thresh = false;    
  
  // loop over readout ranges
  for ( CalDigi::CalXtalReadoutCol::const_iterator it = roCol.begin();
        it !=roCol.end(); it++){

    // get readout range number for both crystal faces
    CalXtalId::AdcRange rangeP = 
      (CalXtalId::AdcRange)it->getRange(CalXtalId::POS); 
    CalXtalId::AdcRange rangeM = 
      (CalXtalId::AdcRange)it->getRange(CalXtalId::NEG); 

    // get adc values 
    int adcP = it->getAdc(CalXtalId::POS);   
    int adcM = it->getAdc(CalXtalId::NEG);   

    float ene;

    // convert adc values into energy
    sc = m_xtalEneTool->calculate(xtalId,
                                  rangeP,rangeM,
                                  adcP,  adcM,
                                  ene, below_thresh);
    if (sc.isFailure()) return sc;

    // create output object
    CalXtalRecData::CalRangeRecData* rangeRec =
      new CalXtalRecData::CalRangeRecData(rangeP,ene,rangeM,ene);

    // add output object to output collection
    recData->addRangeRecData(*rangeRec);
    delete rangeRec;
  }             
  
  return StatusCode::SUCCESS;
}


StatusCode CalXtalRecAlg::computePosition(CalXtalRecData* recData, const CalDigi* digi)

  // Purpose and method:
  //              This function calculates the longitudinal position
  //              for each crystal from light asymmetry.
  //              
  //  Input: CalXtalRecData* recData - pointer to the reconstructed crystal data
  //  Output: the same object, the calculated position is stored using
  //             public function SetPosition()  

{
  MsgStream msg(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;

  // get crystal identification 
  CalXtalId xtalId = digi->getPackedId();

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
  for(int ifield = 0; ifield<fSegment; ifield++)segm11Id.append(segm0Id[ifield]);
  segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
  //get 3D transformation for the last segment of this crystal
  m_detSvc->getTransform3DByID(segm11Id,&transf);
  //get position of the center of the last segment
  Vector vect11 = transf.getTranslation();

  Point p0(0.,0.,0.);           
  // position of the crystal center
  Point pCenter = p0+(vect0+vect11)*0.5; 
  //normalized vector of the crystal direction 
  Vector dirXtal = 0.5*(vect11-vect0)*m_nCsISeg/(m_nCsISeg-1);  

  //-- GET BEST RANGE READOUT --//
  // readout 0 will contain the best ranges in 1-range mode
  // and ALL-RANGE mode
  // theoretically it wont work in 4-range readout mode w/ out
  // range selection, but that mode is rarely used.
  const CalDigi::CalXtalReadoutCol& roCol = digi->getReadoutCol();
  const CalDigi::CalXtalReadout &ro = *roCol.begin();

  int adcP = ro.getAdc(CalXtalId::POS);
  int adcN = ro.getAdc(CalXtalId::NEG);
  CalXtalId::AdcRange rngP = (CalXtalId::AdcRange)ro.getRange(CalXtalId::POS);
  CalXtalId::AdcRange rngN = (CalXtalId::AdcRange)ro.getRange(CalXtalId::NEG);

  float pos;
  sc = m_xtalPosTool->calculate(xtalId,rngP,rngN,adcP,adcN,pos);
  if (sc.isFailure()) return sc;


  // put 1D position info into 3D vector
  // 'pos' is in units of xtal Length, convert to rel units (-1->1)
  pos /= m_CsILength; 
  Point pXtal = pCenter+dirXtal*pos;

  // store calculated position in the reconstructed data
  // for the best readout range
  (recData->getRangeRecData(0))->setPosition(pXtal);

  return StatusCode::SUCCESS;
}
