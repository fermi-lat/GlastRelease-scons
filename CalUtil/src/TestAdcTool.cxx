// Include files
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "CLHEP/Random/RandGauss.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TestAdcTool.h"

// to access an XML containing Digi parameters file
#include "xml/IFile.h"

static ToolFactory<TestAdcTool> s_factory;
const IToolFactory& TestAdcToolFactory = s_factory;

TestAdcTool::TestAdcTool( const std::string& type,
                          const std::string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent) {
  declareInterface<ICalAdcTool>(this);

  declareProperty ("doFluctuations", m_doFluctuations=true);
}

StatusCode TestAdcTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc = StatusCode::SUCCESS;

  double value;
  typedef std::map<int*,std::string> PARAMAP;

  PARAMAP param;
  param[&m_xNum]        = std::string("xNum");
  param[&m_eTowerCal]   = std::string("eTowerCAL");
  param[&m_eLatTowers]  = std::string("eLATTowers");
  param[&m_nCsISeg]     = std::string("nCsISeg");
  param[&m_eXtal]       = std::string("eXtal");
  param[&m_eDiodeMSmall]= std::string("eDiodeMSmall");
  param[&m_eDiodePSmall]= std::string("eDiodePSmall");
  param[&m_eDiodeMLarge]= std::string("eDiodeMLarge");
  param[&m_eDiodePLarge]= std::string("eDiodePLarge");
  param[m_noise]        = std::string("cal.noiseLarge");
  param[m_noise+1]      = std::string("cal.noiseSmall");
  param[m_ePerMeV+1]    = std::string("cal.ePerMeVSmall");
  param[m_ePerMeV]      = std::string("cal.ePerMevLarge");
  param[&m_pedestal]    = std::string("cal.pedestal");
  param[&m_maxAdc]      = std::string("cal.maxAdcValue");

  // now try to find the GlastDevSvc service

  IGlastDetSvc* detSvc;
  sc = service("GlastDetSvc", detSvc);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to get GlastDetSvc " << endreq;
    return sc;
  }

  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){
    if(!detSvc->getNumericConstByName((*it).second, &value)) {
      msglog << MSG::ERROR << " constant " <<(*it).second <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*it).first)=(int)value;
  }

  // extracting double constants from detModel. Store into local cache - member variables.

  typedef std::map<double*,std::string> DPARAMAP;
  DPARAMAP dparam;
  dparam[&m_CsILength]  = std::string("CsILength");
  dparam[&m_lightAtt]   = std::string("cal.lightAtt");
  dparam[&m_thresh]=std::string("cal.zeroSuppressEnergy");


  for(DPARAMAP::iterator dit=dparam.begin(); dit!=dparam.end();dit++){
    if(!detSvc->getNumericConstByName((*dit).second,(*dit).first)) {
      msglog << MSG::ERROR << " constant " <<(*dit).second << " not defined" << endreq;
      return StatusCode::FAILURE;
    }
  }

  // eperMevInDiode originally from CalDigi XML file
  // but i don't want dependency, so i'm hard coding it.  (just a testing tool anyway :)
  m_ePerMeVinDiode = 2.77e5;

  // from CalUtil::LinearConvertAdc
  if(!detSvc->getNumericConstByName("cal.maxResponse0",&m_maxResponse[0])) {
    msglog << MSG::ERROR << " constant " << " cal.maxResponse0 not defined" << endreq;
    return StatusCode::FAILURE;
  }
  if(!detSvc->getNumericConstByName("cal.maxResponse1",&m_maxResponse[1])) {
    msglog << MSG::ERROR << " constant " << " cal.maxResponse1 not defined" << endreq;
    return StatusCode::FAILURE;
  }
  if(!detSvc->getNumericConstByName("cal.maxResponse2",&m_maxResponse[2])) {
    msglog << MSG::ERROR << " constant " << " cal.maxResponse2 not defined" << endreq;
    return StatusCode::FAILURE;
  }
  if(!detSvc->getNumericConstByName("cal.maxResponse3",&m_maxResponse[3])) {
    msglog << MSG::ERROR << " constant " << " cal.maxResponse3 not defined" << endreq;
    return StatusCode::FAILURE;
  }

  for (int face=0; face < 2; face++) {
    for (int range = 0; range < 4; range++) {
      m_gain[face][range] = (m_maxAdc-m_pedestal)/m_maxResponse[range];
    }
  }

  return StatusCode::SUCCESS;
}

/*!
  Method:
  -# loop through hits, make sure they are for appropriate volume
  -# for each hit, add energy to either A) one diode volume (direct diode hit) or B) calculate signal to each face (xtal hit)
  -# sum signal to each diode
  -# calculate random noise and add that to diode signal
  -# calculate adc response for each diode/range adc = diodeEnergy*gain + ped
  -# find best range for each face
*/
StatusCode TestAdcTool::calculate(const idents::CalXtalId &xtalId,
                                  const std::vector<const Event::McIntegratingHit*> &hitList, // list of all mc hits for this xtal & it's diodes.
                                  bool &lacP,
                                  bool &lacN,
                                  idents::CalXtalId::AdcRange &rangeP, // output - best range
                                  idents::CalXtalId::AdcRange &rangeN,  // output - best range
                                  std::vector<int> &adcP,              // output - ADC's for all ranges 0-3
                                  std::vector<int> &adcN              // output - ADC's for all ranges 0-3
                                  ) {

  MsgStream msglog(msgSvc(), name());

  // STOLEN from CalDigiAlg::XtalSignal
  // signal for both xtal faces (POS, NEG)
  std::vector<double> signal(2,0);

  // direct energy depositions in 4 diodes of one xtal; vector contains all 4 diodes
  typedef idents::CalXtalId CalXtalId;
  typedef std::pair<int,int> diodeId;
  std::map<diodeId,double> diodeEnergy;
  // initialize 4 diodes
  diodeEnergy[diodeId(0,0)] = 0;
  diodeEnergy[diodeId(0,1)] = 0;
  diodeEnergy[diodeId(1,0)] = 0;
  diodeEnergy[diodeId(1,1)] = 0;

  // STAGE 1 Clone CalDigiAlg::fillSignalEnergies() ////////////////////////////////////////

  // loop over hits.
  for (std::vector<const Event::McIntegratingHit*>::const_iterator it = hitList.begin();
       it != hitList.end(); it++) {

    // get volume identifier.
    idents::VolumeIdentifier volId = ((idents::VolumeIdentifier)(*it)->volumeID());

    // skip any hits not registered as Cal
    if ((int)volId[fLATObjects] != m_eLatTowers ||
        (int)volId[fTowerObjects] != m_eTowerCal) {
      msglog << MSG::WARNING << "Invalid volume id.  Continuing w/out it but this shouldn't happen." << endreq;
      continue;
    }

    // make sure volumeid matches xtalId
    if ((int)volId[fCALXtal]  != xtalId.getColumn() ||
        (int)volId[fLayer]    != xtalId.getLayer()  ||
        (int)(volId[fTowerY]*m_xNum + volId[fTowerX]) != xtalId.getTower()) {
      msglog << MSG::WARNING << "Invalid volume id.  Continuing w/out it but this shouldn't happen." << endreq;
      continue;
    }

    double ene = (*it)->totalEnergy();
    HepPoint3D mom1 = (*it)->moment1();
    int segm = volId[fSegment];

    // add energy to appropriate volume
    if((int)volId[fCellCmp] == m_eDiodePLarge)
      diodeEnergy[diodeId(CalXtalId::POS,CalXtalId::LARGE)] += ene;

    else if((int)volId[fCellCmp] == m_eDiodeMLarge)
      diodeEnergy[diodeId(CalXtalId::NEG,CalXtalId::LARGE)] += ene;

    else if((int)volId[fCellCmp] == m_eDiodePSmall )
      diodeEnergy[diodeId(CalXtalId::POS,CalXtalId::SMALL)] += ene;

    else if((int)volId[fCellCmp] ==  m_eDiodeMSmall)
      diodeEnergy[diodeId(CalXtalId::NEG,CalXtalId::SMALL)] += ene;

    else if((int)volId[fCellCmp] ==  m_eXtal ){
      // let's define the position of the segment along the crystal
      double relpos = (segm+0.5)/m_nCsISeg;
      // in local reference system x is always oriented along the crystal
      double dpos =  mom1.x();

      relpos += dpos/m_CsILength;

      // take into account light tapering
      std::pair<double,double> signals;
      calculateSignals(relpos,ene,signals);

      // set up a XtalMap to act as key for the map - if there is no entry, add
      // add one, otherwise add signal to existing map element.
      signal[idents::CalXtalId::POS] += signals.first;
      signal[idents::CalXtalId::NEG] += signals.second;
    }
  }

  // STAGE 2 Clone CalDigiAlg::addNoiseToSignals() and ::addNewNoseToSignals() 
  // Purpose and Method:
  // add electronic noise to the diode response and add to the crystal
  // response. The diodeEnergy becomes the readout source.

  // loop through all 4 diodes
  for (int face = 0; face < 2; face++)
    for (int diode = 0; diode < 2; diode++) {
      diodeId diode_id(face,diode);

      // convert energy deposition in a diode to
      // the equivalent energy in a crystal
      diodeEnergy[diode_id] *= m_ePerMeVinDiode/m_ePerMeV[diode];

      // add crystal signal - now diode energy contains
      // the signal at given diode in energy units
      // (equivalent energy deposition at the crystal center)
      diodeEnergy[diode_id] += signal[face];

      // add poissonic fluctuations in the number of electrons in a diode
      if (m_doFluctuations) {
        float numberElectrons = diodeEnergy[diode_id] * m_ePerMeV[diode];

        // approximate Poisson distribution by gaussian for numberElectrons >> 1
        float electronFluctuation = sqrt(numberElectrons) * RandGauss::shoot();

        diodeEnergy[diode_id] += electronFluctuation /m_ePerMeV[diode];
      }

      // add electronic noise
      diodeEnergy[diode_id] += RandGauss::shoot()*m_noise[diode]/m_ePerMeV[diode];
    }

  // STAGE 3 Clone CalDigiAlg::createDigis

  // calculate all ADC responses
  CalXtalId::XtalFace xtalface;
  for (int range = 0; range < 4; range++) {
    int diode = range/2;
    xtalface = CalXtalId::POS;
    diodeId dId(xtalface,diode);
    adcP[range] = (int)(diodeEnergy[dId]*m_gain[xtalface][range]) + m_pedestal;

    xtalface = CalXtalId::NEG;
    dId = diodeId(xtalface,diode);
    adcN[range] = (int)(diodeEnergy[dId]*m_gain[xtalface][range]) + m_pedestal;
    if (adcN[range] > m_maxAdc) adcN[range] = m_maxAdc;
  }

  // Best Range And Lac Tests //

  // best range POS
  xtalface = CalXtalId::POS;
  for (int range = 0; range < 4; range++) {
    int diode = range/2;
    double ene = diodeEnergy[diodeId(xtalface,diode)];
    rangeP = (CalXtalId::AdcRange)range;
    if (ene < m_maxResponse[range]) break;  // break on 1st energy range that is < threshold
  }

  // log-accept POS - hardware only tests large diode for zero suppression
  if (diodeEnergy[diodeId(xtalface,idents::CalXtalId::LARGE)] > m_thresh) 
    lacP = true;
  else lacP = false;

  // best range NEG
  xtalface = CalXtalId::NEG;
  for (int range = 0; range < 4; range++) {
    int diode = range/2;
    double ene = diodeEnergy[diodeId(xtalface,diode)];
    rangeN = (CalXtalId::AdcRange)range;
    if (ene < m_maxResponse[range]) break;  // break on 1st energy range that is < threshold
  }

  // log-accept NEG - hardware only tests large diode fro zero suppression
  lacN = false;
  if (diodeEnergy[diodeId(xtalface,idents::CalXtalId::LARGE)] > m_thresh) 
    lacN = true;
  else lacN = false;

  msglog << MSG::DEBUG  << "CalAdcTool id " << xtalId
         << " s0=" << signal[CalXtalId::POS]
         << " s1=" << signal[CalXtalId::NEG]
         << " d0=" << diodeEnergy[diodeId(0,0)]
         << " d1=" << diodeEnergy[diodeId(0,1)]
         << " d2=" << diodeEnergy[diodeId(1,0)]
         << " d3=" << diodeEnergy[diodeId(1,1)]
         << endreq;

  return StatusCode::SUCCESS;
}

void TestAdcTool::calculateSignals(double relativePosition,
                                   double depositedEnergy,
                                   std::pair<double, double> &signals) {

  // Purpose and Method: calculate light taper for the two crystal ends .
  // Inputs: crystal id, deposited energy and relative position in the crystal
  // Outputs: energy seen at the two crystal ends
  // stolen from CalDigi by Z. Fewtrell 8/04

  double norm = 0.5+0.5*m_lightAtt; // light tapering in the center of crystal (relpos=0.5)

  double s2 = depositedEnergy*(1-relativePosition*(1-m_lightAtt))/norm;
  double s1 = depositedEnergy*(1-(1-relativePosition)*(1-m_lightAtt))/norm;

  signals = std::make_pair(s1,s2);
};
