//    $Header$

/** @file
    @author Z.Fewtrell
*/

// Include files

// LOCAL
#include "INeighborXtalkTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "Event/Digi/CalDigi.h"
#include "idents/CalXtalId.h"
#include "facilities/Util.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"


// STD
#include <string>
#include <map>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace CalUtil;
using namespace Event;
using namespace idents;
using namespace facilities;

/// used to indicate empty channel
static const short INVALID_ADC = -5000;


/** \brief Simple implementation of INeighborXtalkTool
    \note currently reads xtalk curves from txt file & not from calib db.
    \note currently only deals w/ cross-talk listed as he_diode->he_diode

    jobOptions:
    - txtFile - (default="$(CALXTALRESPONSEDATAPATH)/Xtalk/CU06_Neighbor_xtalk.txt") input xtalk coefficients
    - CalCalibSvc - (default="CalCalibSvc") - source for Cal Calibrations
*/
class NeighborXtalkTool : public AlgTool, 
                          virtual public INeighborXtalkTool {
public:
  /// default ctor, declares jobOptions.
  NeighborXtalkTool( const string& type, 
                     const string& name, 
                     const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize() {return StatusCode::SUCCESS;}


  StatusCode calcXtalkCIDAC(CalUtil::DiodeIdx diodeIdx, float &xtalkCIDAC) const ;
  
  StatusCode calcXtalkMeV(CalUtil::DiodeIdx diodeIdx, float &xtalkMev) const;

  StatusCode buildSignalMap(const Event::CalDigiCol &digiCol);

private:
  /// load in xtalk table from txt ifle
  StatusCode loadXtalkTable(const string &infile);

  /// represents single xtalk curve
  struct XtalkEntry {
    float x_intercept;
    float slope;
  };

  /// type for associating single adc channel w/ xtalk curve
  typedef map<FaceIdx, XtalkEntry> ChannelSplineMap;

  /// type for associating (dest_channel,src_channel) tuple w/ xtalk curve
  typedef map<FaceIdx, ChannelSplineMap > XtalkMap;

  void insertXtalkCurve(FaceIdx srcIdx,
                        FaceIdx destIdx,
                        const XtalkEntry &xtalk);

  /// store xtalk curves associated w/ channels.
  XtalkMap m_xtalkMap;

  /// read in xtalk table from this file.
  StringProperty m_txtFilename;

  /// \brief store output signals for HE diode on each cal face
  ///  signal stored in he_dac units
  CalVec<FaceIdx, float> m_signalMap;


  /// name of CalCalibSvc to use for calib constants.
  StringProperty m_calCalibSvcName;                         
  /// pointer to CalCalibSvc object.
  ICalCalibSvc *m_calCalibSvc;  

  /// evaluate xtalk from single channel
  /// xtalk in dac units, retuns 0 if unapplicable
  float evalSingleChannelXtalk(float srcDac, const XtalkEntry &xtalk) const;

};

//static ToolFactory<NeighborXtalkTool> s_factory;
//const IToolFactory& NeighborXtalkToolFactory = s_factory;
DECLARE_TOOL_FACTORY(NeighborXtalkTool);

NeighborXtalkTool::NeighborXtalkTool( const string& type, 
                                      const string& name, 
                                      const IInterface* parent) :
  AlgTool(type,name,parent),
  m_signalMap(FaceIdx::N_VALS, INVALID_ADC),
  m_calCalibSvc(0)
{
  declareInterface<INeighborXtalkTool>(this);

  declareProperty("txtFile", m_txtFilename="$(CALXTALRESPONSEDATAPATH)//Xtalk/CU06_Neighbor_xtalk.txt");
  declareProperty("CalCalibSvc", m_calCalibSvcName="CalCalibSvc");
  
}

StatusCode NeighborXtalkTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }
  
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName  << endreq;
    return sc;
  }
  
  //-- read xtalk curves table from disk --//
  if ((sc = loadXtalkTable(m_txtFilename)).isFailure())
    return sc;

  return StatusCode::SUCCESS;
}

StatusCode NeighborXtalkTool::calcXtalkCIDAC(CalUtil::DiodeIdx diodeIdx, float &xtalkCIDAC) const {
  // start w/ 0 xtalk in case of early function exit.
  xtalkCIDAC = 0;
  /// currently only simulating crosstalk w/ HE diode destination channe,l
  if (diodeIdx.getDiode() != SM_DIODE)
    return StatusCode::SUCCESS;

  XtalkMap::const_iterator xtalkIt = m_xtalkMap.find(diodeIdx.getFaceIdx());

  // return 0 if there are no neighboring source channels w/ registered xtalk
  if (xtalkIt == m_xtalkMap.end())
    return StatusCode::SUCCESS;

  // find all source curves for given destination channel
  for (ChannelSplineMap::const_iterator chanIt =
         xtalkIt->second.begin();
       chanIt != xtalkIt->second.end();
       chanIt++) {
    FaceIdx srcIdx(chanIt->first);
    float srcDac = m_signalMap[srcIdx];

    if (srcDac <= 0)
      continue;

    xtalkCIDAC += evalSingleChannelXtalk(srcDac, chanIt->second);
  }

  return StatusCode::SUCCESS;
}



StatusCode NeighborXtalkTool::loadXtalkTable(const string &filename) {
  // open file
  string filenameCopy(filename);
  Util::expandEnvVar(&filenameCopy);
  MsgStream msglog(msgSvc(), name());

  msglog << MSG::INFO << "Loading xtalk tables from: " << filename << endreq;

  ifstream infile(filenameCopy.c_str());
  if (!infile.is_open()) {
    msglog << MSG::ERROR << "can't find xtalk table file: " << filename << endreq;

    return StatusCode::FAILURE;
  }

  // loop through each line in file
  string line;
  while (infile.good()) {
    getline(infile, line);
    if (infile.fail()) break; // bad get

    // check for comments
    if (line[0] == ';')
      continue;

    istringstream istrm(line);

    // read values from line
    unsigned srcDiodeIdx, destDiodeIdx;
    float x_int, slope;
    istrm >> destDiodeIdx >> srcDiodeIdx >> x_int >> slope;

    DiodeIdx srcIdx, destIdx;
    srcIdx.setVal(srcDiodeIdx);
    destIdx.setVal(destDiodeIdx);

    /// currently only simulating crosstalk w/ HE diode -> HE diode
    if (srcIdx.getDiode() != SM_DIODE)
      continue;
    if (destIdx.getDiode() != SM_DIODE)
      continue;

    XtalkEntry xtalk = {x_int, slope};
        
    insertXtalkCurve(srcIdx.getFaceIdx(),
                     destIdx.getFaceIdx(),
                     xtalk);
  }


  return StatusCode::SUCCESS;
}


void NeighborXtalkTool::insertXtalkCurve(FaceIdx srcIdx,
                                         FaceIdx destIdx,
                                         const XtalkEntry &xtalk) {
  // find all cross talk entries for given 'destination' channel
  XtalkMap::iterator xtalkIt = m_xtalkMap.find(destIdx);


  // create new destination map if needed
  if (xtalkIt == m_xtalkMap.end())
    xtalkIt = m_xtalkMap.insert(XtalkMap::value_type(destIdx, ChannelSplineMap())).first;

  // find curve for given source, destination pair.
  ChannelSplineMap::iterator chanIt =
    xtalkIt->second.find(srcIdx);

  // create new spline curve if needed
  if (chanIt == xtalkIt->second.end())
    chanIt = xtalkIt->second.insert(ChannelSplineMap::value_type(srcIdx, xtalk)).first;
}

StatusCode NeighborXtalkTool::buildSignalMap(const Event::CalDigiCol &digiCol) {
  // re-initialize signalMap
  fill(m_signalMap.begin(), m_signalMap.end(), INVALID_ADC);

  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator digiIter = digiCol.begin(); 
       digiIter != digiCol.end(); digiIter++) {

    // if there is no digi data, then move on w/ out creating
    // recon TDS data for this xtal
    if ((*digiIter)->getReadoutCol().size() < 1) continue;
    
    CalXtalId xtalId((*digiIter)->getPackedId());
    XtalIdx xtalIdx(xtalId);

    // currently always using 1st readout
    CalDigi::CalXtalReadoutCol::const_iterator ro = 
      (*digiIter)->getReadoutCol().begin();

    for (FaceNum face; face.isValid(); face++) {
      RngNum rng(ro->getRange(face));
      if (rng.getDiode() == LRG_DIODE)
        continue;
      const float adc(ro->getAdc(face));

      const FaceIdx faceIdx(xtalIdx,face);
      const RngIdx rngIdx(faceIdx,rng);

      // pedestal subtract
      float ped;
      StatusCode sc = m_calCalibSvc->getPed(rngIdx,ped);
      if (sc.isFailure()) {
        MsgStream msglog(msgSvc(), name());
        msglog << MSG::ERROR << "can't find cal ped for given digi channel: " << rngIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
      const float adcPed = adc - ped;
          
      // evaluate CIDAC signal
      float cidac = 0;
      sc = m_calCalibSvc->evalCIDAC(rngIdx, adcPed, cidac);
      if (sc.isFailure()) return sc;

      m_signalMap[faceIdx] = cidac;
    }
  }
  
  return StatusCode::SUCCESS;
}


float NeighborXtalkTool::evalSingleChannelXtalk(const float srcDac, const XtalkEntry &xtalk)  const {
  if (srcDac <= xtalk.x_intercept)
    return 0;
  
  // currently simple linear model
  return (srcDac - xtalk.x_intercept)*xtalk.slope;
}

StatusCode NeighborXtalkTool::calcXtalkMeV(const CalUtil::DiodeIdx diodeIdx, float &xtalkMev) const {
  xtalkMev = 0;
  StatusCode sc;
        
  float xtalkCIDAC;
  sc = calcXtalkCIDAC(diodeIdx, xtalkCIDAC);
  if (sc.isFailure()) return sc;

  // quick exit.
  if (xtalkCIDAC == 0)
    return StatusCode::SUCCESS;

  float mpdDiode;
  sc = m_calCalibSvc->getMPDDiode(diodeIdx, mpdDiode);
  if (sc.isFailure()) return sc;

  /// mev = (mev/cidac)*cidac
  xtalkMev = xtalkCIDAC*mpdDiode;

  return StatusCode::SUCCESS;
}
