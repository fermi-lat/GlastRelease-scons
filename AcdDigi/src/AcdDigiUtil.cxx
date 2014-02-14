#define AcdDigi_AcdDigiUtil_CPP 

// File and Version Information:
// $Header$
// Description
// Some utility methods helpful for performing the ACD digitization.


#include "AcdDigiUtil.h"

#include "GaudiKernel/MsgStream.h"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdCoherentNoise.h"
#include "CalibData/Acd/AcdRibbon.h"

#include "AcdUtil/AcdCalibFuncs.h"
#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdGeomMap.h"
#include "AcdUtil/IAcdCalibSvc.h"

#include "Event/MonteCarlo/McPositionHit.h"

// for min and floor functions
#include <algorithm>
#include <cmath>

#include <map>
#include <utility>

// to access an XML containing Digi parameters file
#include "xmlBase/IFile.h"
#include "facilities/Util.h"

AcdSimCalibData::AcdSimCalibData()
  :m_ped(0),
   m_gain(0),
   m_highRange(0),
   m_range(0),  
   m_veto(0),
   m_cno(0),
   m_coherentNoise(0),
   m_ribbon(0),
   m_pe_per_mip(0.),
   m_mip_per_MeV(0.),
   m_threshold_pha(0),    
   m_veto_threshold_mips(0.),
   m_xover_mips(0.),
   m_cno_threshold_mips(0.),
   m_pedestal_highRange(0.),
   //Change made by D. Green to incorporate zero suppression as a JO
   m_pedestal(0.)
{;}

StatusCode AcdSimCalibData::latchPePerMeV( ) {
  m_pe_per_MeV = m_pe_per_mip * m_mip_per_MeV;
  return StatusCode::SUCCESS;
}

StatusCode AcdSimCalibData::latchPhaThreshold(double countsAbovePed) {
  if  ( m_ped == 0 ) return StatusCode::FAILURE;
  m_threshold_pha = (unsigned short) ( m_ped->getMean() + countsAbovePed );
  return StatusCode::SUCCESS;
}

//Change made by D. Green to incorporate zero suppression as a JO
StatusCode AcdSimCalibData::latchPedestal() {
  if  ( m_ped == 0 ) return StatusCode::FAILURE;
  m_pedestal = (unsigned short) ( m_ped->getMean());
  return StatusCode::SUCCESS;
}

StatusCode AcdSimCalibData::latchVetoThreshold() {
  if ( m_ped == 0 || m_gain == 0 || m_veto == 0 ) return StatusCode::FAILURE;
  double pedestal = m_ped->getMean();
  double mipPeak = m_gain->getPeak();
  unsigned short veto_pha = (unsigned short)(m_veto->getVeto());
  return AcdCalib::mipEquivalent_lowRange(veto_pha,pedestal,mipPeak,m_veto_threshold_mips);
}

StatusCode AcdSimCalibData::latchCnoThreshold() {
  if ( m_cno == 0 || m_highRange == 0 ) return StatusCode::FAILURE;
  double pedestal = m_highRange->getPedestal();
  double slope = m_highRange->getSlope();
  double saturation = m_highRange->getSaturation();  
  unsigned short cno_pha = (unsigned short)(m_cno->getCno());
  return AcdCalib::mipEquivalent_highRange(cno_pha,pedestal,slope,saturation,m_cno_threshold_mips);
}

StatusCode AcdSimCalibData::latchXOverMips() {
  if ( m_ped == 0 || m_gain == 0 || m_range == 0 ) return StatusCode::FAILURE;
  double pedestal = m_ped->getMean();
  double mipPeak = m_gain->getPeak();
  unsigned short xover_pha = (unsigned short)(m_range->getLowMax());
  return AcdCalib::mipEquivalent_lowRange(xover_pha,pedestal,mipPeak,m_xover_mips); 
}


double AcdDigiUtil::simulateDynodeChain(double pe_meanValue)  {
 
  if ( pe_meanValue <= 0 ) return 0.;
  
  unsigned poisson_steps = 4;  //Number of "dynodes" considered
  double gain=5.;         //Nominal gain in each dynode

  double norm(1.);
  double val = pe_meanValue;
  for ( unsigned i(0); i < poisson_steps; i++ ) {
    val = shootPoisson(val);
    if ( val <= 0. ) return val;
    val *= gain;
    norm *= gain;
    // no need to simulate huge signals
    if ( val > 1000. ) break;
  }
  
  val /= norm;
  return val;
}

double AcdDigiUtil::shootPoisson(double pmtPhotoElectrons) {
    // Pupose and Method:  Returns a value from a Poisson distribution,
    //   using the input number of photoelectrons as the mean
    // Input:
    // Output: 
    
    return CLHEP::RandPoisson::shoot(pmtPhotoElectrons);
}

double AcdDigiUtil::shootGaussian(double std_dev) {
    // Purpose and Method:  Returns a value from a gaussian distribution, with mean
    //   zero and standard deviation determined by the one input parameter.
    // Input:  Standard Deviation
    // Output:  A value obtained from the the Gaussian distribution.
    
    return CLHEP::RandGauss::shoot(0.0, std_dev);
}

bool AcdDigiUtil::compareVolIds(const idents::VolumeIdentifier& tileId, 
                                const idents::VolumeIdentifier& screwVolId) {

  ///Purpose and Method:  Compare the volumeIdentfiers of a tile and a 
  /// AcdScrewSq to see if the AcdScrewSq occurs within the given tile
  
  if ((tileId[0] != 1) && (screwVolId[0] != 1) ) return false;
  unsigned int i;
  // compare the entries one by one to see if they are equal
  for (i = 0; i<5; i++) {
    if(tileId[i] != screwVolId[i]) return false;
  }
  return true;
}

bool AcdDigiUtil::getRibbonLengthAndBin(const Event::McPositionHit* hit,
					IAcdGeometrySvc& geomSvc,
					MsgStream& log,
					double& ribbonLength, int& ribbonBin) {

  idents::VolumeIdentifier volId = hit->volumeID();
  idents::AcdId ribId(volId);  

  const AcdRibbonDim* ribbonDim = geomSvc.geomMap().getRibbon(ribId,geomSvc);
  if ( ribbonDim == 0 ) {
    log << MSG::ERROR << "no ribbon dim " << volId.name() << std::endl;
    return false;
  }

  // In local coordinates the box should be centered at (0,0,0)
  const HepPoint3D geant_x0 = hit->entryPoint();

  return ribbonDim->getRibbonLengthAndBin(volId,geant_x0,ribbonLength,ribbonBin);
}
			  

AcdDigiUtil::AcdDigiUtil()
  :m_ifile(0){;}


AcdDigiUtil::~AcdDigiUtil() {
  if (m_ifile) delete m_ifile;
}


StatusCode AcdDigiUtil::initialize(AcdUtil::IAcdCalibSvc& calibSvc, IAcdGeometrySvc& geomSvc, 
				   const std::string &xmlFileName) {
  // Purpose and Method:  Read in the parameters from the XML file using IFile
  
    
  // latch the services
  m_calibSvc = &calibSvc;
  m_geomSvc = &geomSvc;

  m_ifile = new xmlBase::IFile(xmlFileName.c_str());
  
  // pe per MIP
  m_mean_pe_per_mip = m_ifile->getDouble("global_constants", "mean_pe_per_mip",18.0);
  m_mean_pe_per_mip_ribbon = m_ifile->getDouble("global_constants", "mean_pe_per_mip_ribbon",4.5);
    
  // mip per MeV
  m_mip_per_MeV = m_ifile->getDouble("global_constants", "mip_per_MeV",0.52631);
  m_mip_per_MeV_ribbon = m_ifile->getDouble("global_constants", "mip_per_MeV_ribbon",0.52631);

  // counts above pedestal to set PHA threshold
  m_counts_above_pedestal = m_ifile->getDouble("global_constants", "counts_above_pedestal", 15);

  m_veto_mips = m_ifile->getDouble("global_constants", "veto_mips",-1);
  m_cno_mips = m_ifile->getDouble("global_constants", "cno_mips",-1);
  
  // noise levels
  m_veto_width_mips = m_ifile->getDouble("global_constants", "noise_std_dev_veto", 0.02);
  m_cno_width_mips = m_ifile->getDouble("global_constants", "noise_std_dev_cno", 0.02);

  // Edge Effects
  m_max_edge_dist = m_ifile->getDouble("edge_effects", "max_edge_dist", 20.0);
  m_edge_slope = m_ifile->getDouble("edge_effects", "edge_slope", 0.01);  
  m_edge_intercept = m_ifile->getDouble("edge_effects", "edge_intercept", 0.8);

  return StatusCode::SUCCESS;

}

StatusCode AcdDigiUtil::photoElectronsFromEnergy(const idents::AcdId& id, double energy,
						 double& pe_pmtA, double& pe_pmtB) {
  
  StatusCode sc = StatusCode::SUCCESS;

  AcdSimCalibData* pmtACalib(0);
  AcdSimCalibData* pmtBCalib(0);
  sc = getCalibData(id,pmtACalib,pmtBCalib);
  if ( sc.isFailure() ) return sc;

  // check the map
  double pePerMeV_A = pmtACalib->pe_per_MeV();
  double pePerMeV_B = pmtBCalib->pe_per_MeV();

  // OK, have the numbers, now just combine them
  sc = AcdCalib::photoElectronsFromEnergy(energy,pePerMeV_A,pe_pmtA);
  if ( sc.isFailure() ) return sc;

  sc = AcdCalib::photoElectronsFromEnergy(energy,pePerMeV_B,pe_pmtB);
  if ( sc.isFailure() ) return sc;

  return sc;
}

StatusCode AcdDigiUtil::photoElectronsFromEnergy_tile(const Event::McPositionHit *hit, bool edgeEffect, 
						      MsgStream& log,
						      double& pe_pmtA, double& pe_pmtB) {

  StatusCode sc = StatusCode::SUCCESS;

  double energy = hit->depositedEnergy();
  
  // check the edge effect
  if ( edgeEffect ) {
    sc = tileEdgeEffect(hit,log,energy);
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "Couldn't apply edge effect " << hit->volumeID().name() << ' ' << energy << endreq;
      return sc;
    }
  }

  double energy_A = energy;
  double energy_B = energy;
  
  // Added in Feb. 2014 by EAC to allow for single hits coming from overlays
  if ( ( hit->getPackedFlags() & 0xF0000000 ) == (unsigned int)Event::AcdDigi::DIGI_OVERLAY ) { 
    bool acceptPMT_A = (hit->getPackedFlags() & 0x0000001) != 0;
    bool acceptPMT_B = (hit->getPackedFlags() & 0x0001000) != 0;
    if ( acceptPMT_A && ( ! acceptPMT_B ) ) {
      energy_A += energy_B;
      energy_B = 0.;
    }
    if ( acceptPMT_B && ( ! acceptPMT_A ) ) {
      energy_B += energy_A;
      energy_A = 0.;
    }
  }

  idents::VolumeIdentifier volId = hit->volumeID();
  idents::AcdId tileId(volId);

  AcdSimCalibData* pmtACalib(0);
  AcdSimCalibData* pmtBCalib(0);
  sc = getCalibData(tileId,pmtACalib,pmtBCalib);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "Couldn't get calib data " << tileId.id() << ' ' << energy << endreq;
    return sc;
  }

  // check the map
  double pePerMeV_A = pmtACalib->pe_per_MeV();
  double pePerMeV_B = pmtBCalib->pe_per_MeV();

  // OK, have the numbers, now just combine them
  sc = AcdCalib::photoElectronsFromEnergy(energy_A,pePerMeV_A,pe_pmtA);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "photoElectronsFromEnergy failed for PMT A " << tileId.id() << ' ' << energy << endreq;
    return sc;
  }

  sc = AcdCalib::photoElectronsFromEnergy(energy_B,pePerMeV_B,pe_pmtB);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "photoElectronsFromEnergy failed for PMT B " << tileId.id() << ' ' << energy << endreq;
    return sc;
  }

  return sc;
}

StatusCode AcdDigiUtil::photoElectronsFromEnergy_ribbon(const Event::McPositionHit *hit,
							MsgStream& log,
							double& pe_pmtA, double& pe_pmtB) {
  
  StatusCode sc = StatusCode::SUCCESS;

  double energy = hit->depositedEnergy();
  idents::VolumeIdentifier volId = hit->volumeID();
  idents::AcdId ribId(volId);

  double ribbonLength(0);
  int ribbonBin(-1);
  
  AcdSimCalibData* pmtACalib(0);
  AcdSimCalibData* pmtBCalib(0);
  sc = getCalibData(ribId,pmtACalib,pmtBCalib);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "Couldn't get calib data " << ribId.id() << ' ' << energy << endreq;
    return sc;
  }


  bool ok = getRibbonLengthAndBin(hit,*m_geomSvc,log,ribbonLength,ribbonBin);
  if ( !ok ) {
    log << MSG::ERROR << "getRibbonLengthAndBin failed " << ribId.id() << ' ' << energy << endreq;
    return StatusCode::FAILURE;
  }

  // check the map
  double pePerMeV_A = pmtACalib->pe_per_MeV();
  double pePerMeV_B = pmtBCalib->pe_per_MeV();

  // fold in ribbon attenuation factors
  if ( pmtACalib->ribbon_calib() == 0 ||
       pmtBCalib->ribbon_calib() == 0 ) {
    log << MSG::ERROR << "Missing a ribbon calibration " << ribId.id() << ' ' << energy << endreq;
    return StatusCode::FAILURE;
  }

  
  double factorA(1.);
  double factorB(1.);
  switch ( ribbonBin ) {
  case 0:
  case 1:
  case 2:
    factorA = pmtACalib->ribbon_calib()->operator[](ribbonBin);
    factorB = pmtBCalib->ribbon_calib()->operator[](ribbonBin);
    break;
  case 4:
  case 5:
  case 6:
    factorA = pmtACalib->ribbon_calib()->operator[](ribbonBin-1);
    factorB = pmtBCalib->ribbon_calib()->operator[](ribbonBin-1);
    break;
  }
  pePerMeV_A *= factorA;
  pePerMeV_B *= factorB;  

  // OK, have the numbers, now just combine them
  sc = AcdCalib::photoElectronsFromEnergy(energy,pePerMeV_A,pe_pmtA);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "photoElectronsFromEnergy failed for PMT A " << ribId.id() << ' ' << energy << endreq;
    return sc;
  }

  sc = AcdCalib::photoElectronsFromEnergy(energy,pePerMeV_B,pe_pmtB);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "photoElectronsFromEnergy failed for PMT B " << ribId.id() << ' ' << energy << endreq;
    return sc;
  }

  log << MSG::DEBUG << "Ribbon: " << ribId.id() << ' ' << volId.name() << ' ' << ribbonBin 
      << "  E: " << energy 
      << "  pe/MeV " << pePerMeV_A << ',' << pePerMeV_B
      << "  pe: " << pe_pmtA << ',' << pe_pmtB << endreq;  

  return sc;
}


/// calulates the light yield expressed in mip equivalent from the number of observed PE 
StatusCode AcdDigiUtil::mipEquivalentLightYeild(const idents::AcdId& id, double pe_pmtA, double pe_pmtB,
						MsgStream& log,
						double& mipEquivA, double& mipEquivB) {
    
  double pePerMip_A(0.);
  double pePerMip_B(0.);

  AcdSimCalibData* pmtACalib(0);
  AcdSimCalibData* pmtBCalib(0);

  StatusCode sc = getCalibData(id,pmtACalib,pmtBCalib);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "Couldn't get calib data " << id.id() << endreq;
    return sc;
  }
  pePerMip_A = pmtACalib->pe_per_mip();
  pePerMip_B = pmtBCalib->pe_per_mip();
  
  sc = AcdCalib::lightYeildMipEquivalent(pe_pmtA,pePerMip_A,mipEquivA);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "lightYeildMipEquivalent failed for PMT A " << id.id() << endreq;
    return sc;
  }

  sc = AcdCalib::lightYeildMipEquivalent(pe_pmtB,pePerMip_B,mipEquivB);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "lightYeildMipEquivalent failed for PMT B " << id.id() << endreq;
    return sc;
  }

  return sc;

}

/// get the PHA counts from the mip equivalent light yield
StatusCode AcdDigiUtil::phaCounts(const idents::AcdId& id, const double mipEquiv[2], bool applyNoise, 
				  MsgStream& log,	       
				 Event::AcdDigi::Range range[2], unsigned short pha[2]) {
  
  AcdSimCalibData* pmtACalib(0);
  AcdSimCalibData* pmtBCalib(0);

  StatusCode sc = getCalibData(id,pmtACalib,pmtBCalib);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "Couldn't get calib data " << id.id() << endreq;
    return sc;
  }

  // Kill 1% of PHA values dead, for testing Ninja stuff
  //if ( CLHEP::RandFlat::shoot(1.) < 0.01 ) {
  //  range[0] = range[1] = Event::AcdDigi::LOW;
  //  pha[0] = pha[1] = 0;
  //  return sc;
  //}

  for ( unsigned i(0); i < 2; i++ ) {
    const AcdSimCalibData* calibData = i == 0 ? pmtACalib : pmtBCalib;

    double xOverMips = calibData->xover_mips();
    range[i] = mipEquiv[i] > xOverMips ? Event::AcdDigi::HIGH : Event::AcdDigi::LOW;
  
    if ( range[i] == Event::AcdDigi::LOW ) { 
      double pedestal = calibData->ped_calib()->getMean();
      double mipPeak = calibData->gain_calib()->getPeak();
      sc = AcdCalib::PHA_lowRange(mipEquiv[i],pedestal,mipPeak,pha[i]);
      if ( applyNoise ) {
	short noise = (short)shootGaussian( calibData->ped_calib()->getWidth() );
	pha[i] = ( -1 * noise < pha[i] ) ? pha[i] + noise : 0;
      }
      if ( sc.isFailure() ) {
	log << MSG::ERROR << "Low range PHA conversion failed " << id.id() << endreq;
	return sc;
      }
    } else {
      double slope = calibData->highRange_calib()->getSlope();
      double pedestal = calibData->highRange_calib()->getPedestal();
      double saturation = calibData->highRange_calib()->getSaturation();
      sc = AcdCalib::PHA_highRange(mipEquiv[i],pedestal,slope,saturation,pha[i]);
      if ( sc.isFailure() ) {
	log << MSG::ERROR << "High range PHA conversion failed " << id.id() << endreq;
	return sc;
      }
    }
  }
  return StatusCode::SUCCESS;
}


/// get the PHA counts from the mip equivalent light yield
StatusCode AcdDigiUtil::applyCoherentNoiseToPha(const idents::AcdId& id, unsigned int deltaGemEventTime, MsgStream& log,
						unsigned short pha[2]) {

  if ( deltaGemEventTime > 2500 ) return StatusCode::SUCCESS;

  AcdSimCalibData* pmtACalib(0);
  AcdSimCalibData* pmtBCalib(0);

  StatusCode sc = getCalibData(id,pmtACalib,pmtBCalib);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "Couldn't get calib data " << id.id() << endreq;
    return sc;
  }

  for ( unsigned i(0); i < 2; i++ ) {
    const AcdSimCalibData* calibData = i == 0 ? pmtACalib : pmtBCalib;
    const CalibData::AcdCoherentNoise* cNoise = calibData->coherentNoise_calib();
    double deltaPed(0.);
    sc = AcdCalib::coherentNoise(deltaGemEventTime,
				 cNoise->getAmplitude(),cNoise->getDecay(),cNoise->getFrequency(),cNoise->getPhase(),
				 deltaPed); 
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "Failed to apply coherent noise " << id.id() << endreq;
      return sc;
    }
    if ( -1 * deltaPed > pha[i] ) {
      pha[i] = 0;
    } else {
      pha[i] += (short unsigned int)(deltaPed);
    }
  }
  return StatusCode::SUCCESS;
}

/// Adjusts the deposited energy recorded in an ACD volume 
/// based on the location of the hit
StatusCode AcdDigiUtil::tileEdgeEffect(const Event::McPositionHit *hit, MsgStream& log, double& energy) {
      

  idents::VolumeIdentifier volId = hit->volumeID();
  idents::AcdId tileId(volId);
  
  // In local coordinates the box should be centered at (0,0,0)
  const HepPoint3D geant_x0 = hit->entryPoint();

  const AcdTileDim* tileDim = m_geomSvc->geomMap().getTile(tileId,*m_geomSvc);
  if ( tileDim == 0 ) { 
    log << MSG::ERROR << "Couldn't get tile geometry " << volId.name() << endreq;
    return StatusCode::FAILURE;
  }
  
  AcdFrameUtil::AcdReferenceFrame refFrame = m_geomSvc->getReferenceFrame(volId);
  const HepGeom::Transform3D& rotToLoca = AcdFrameUtil::getRotationToLocal(refFrame);
  const HepPoint3D local_x0 = rotToLoca*geant_x0;
  
  double activeX(0.);
  double activeY(0.);
  
  const AcdTileSection* sect = tileDim->getSection(volId[6]);
  sect->activeDistance(local_x0, volId[6], activeX, activeY );
  double dist = activeX > activeY ? activeY : activeX;
  
    // Apply edge correction if within m_max_edge_dist (mm) of the edge
  if (dist < m_max_edge_dist) {
    energy = (m_edge_slope*dist + m_edge_intercept) * hit->depositedEnergy();   
  } else {
    energy = hit->depositedEnergy();
  }

  return StatusCode::SUCCESS;
}


/// Checks all the various thresholds
///Change made by D. Green to incorporate zero suppression as a JO
StatusCode AcdDigiUtil::checkThresholds(const idents::AcdId& id, const double mipEquiv[2],
					const unsigned short phaArr[2], const Event::AcdDigi::Range rangeArr[2], bool applyNoise, 
					MsgStream& log,
					bool& makeDigi, 
					bool phaThreshArr[2], bool vetoArr[2],  bool highArr[2],
					const double phaZeroThreshold) {
  
  AcdSimCalibData* pmtACalib(0);
  AcdSimCalibData* pmtBCalib(0);

  StatusCode sc = getCalibData(id,pmtACalib,pmtBCalib);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "Couldn't get calib data " << id.id() << endreq;
    return sc;
  }

  makeDigi = false;  
  for ( unsigned i(0); i < 2; i++ ) {
    const AcdSimCalibData* calibData = i == 0 ? pmtACalib : pmtBCalib;
  
    phaThreshArr[i] = rangeArr[i] == Event::AcdDigi::HIGH ? true : (phaArr[i] >= (calibData->pedestal_pha() + phaZeroThreshold));
    double mipVeto = mipEquiv[i];
    double cnoVeto = mipEquiv[i];
    if ( applyNoise ) { 
      mipVeto += shootGaussian( calibData->veto_width_mips() );
      cnoVeto += shootGaussian( calibData->cno_width_mips() );
    }
    vetoArr[i] = mipVeto >= calibData->veto_threshold_mips();
    highArr[i] = cnoVeto >= calibData->cno_threshold_mips();
    makeDigi |= phaThreshArr[i] || vetoArr[i] || highArr[i];
  }

  return StatusCode::SUCCESS;
}			     
  



StatusCode AcdDigiUtil::getCalibData(const idents::AcdId& id, AcdSimCalibData*& pmtACalib, AcdSimCalibData*& pmtBCalib) {

  std::map< unsigned int, std::pair<AcdSimCalibData*,AcdSimCalibData*> >::const_iterator itrFind = m_calibMap.find(id);  
  bool upToDate(false);
  if ( itrFind != m_calibMap.end() ) {
    // get it
    pmtACalib = itrFind->second.first;
    pmtBCalib = itrFind->second.second;
        // IN simulation there should really only be one version of the calibrations
    upToDate = true;
  } else {
    pmtACalib = new AcdSimCalibData;
    pmtBCalib = new AcdSimCalibData;
    // put the pmt calib in the map
    m_calibMap[id] = std::make_pair(pmtACalib,pmtBCalib);
  }    

  StatusCode sc = StatusCode::SUCCESS;
  // possible reflush the data, though this should only happen the first time we get the data
  if ( ! upToDate ) {
    sc = fetchCalibData(id,Event::AcdDigi::A,pmtACalib);
    if ( sc.isFailure() ) return sc;
    sc = fetchCalibData(id,Event::AcdDigi::B,pmtBCalib);
  }

  return sc;
}


StatusCode AcdDigiUtil::fetchCalibData(const idents::AcdId& id, Event::AcdDigi::PmtId pmt, AcdSimCalibData*& pmtCalib ) {

  StatusCode sc = StatusCode::SUCCESS;

  CalibData::AcdPed* ped(0);
  CalibData::AcdGain* gain(0);
  CalibData::AcdHighRange* highRange(0);
  CalibData::AcdRange* range(0);  
  CalibData::AcdVeto* veto(0);
  CalibData::AcdCno* cno(0);
  CalibData::AcdCoherentNoise* coherentNoise(0);  
  CalibData::AcdRibbon* ribbon(0);  
  

  // Check the XML file
  std::string idStr;
  facilities::Util::itoa(id.id(), idStr);
  std::string pmtIdStr = idStr;
  if ( m_ifile->contains("meanPePerMip", pmtIdStr.c_str())){
    std::vector<double> pePerMipVec = m_ifile->getDoubleVector("meanPePerMip", pmtIdStr.c_str());
    pmtCalib->setPe_per_mip( pePerMipVec[pmt] );
  } else {
    if ( id.tile() ) {
      pmtCalib->setPe_per_mip( m_mean_pe_per_mip );
    } else if ( id.ribbon() ) {
      pmtCalib->setPe_per_mip( m_mean_pe_per_mip_ribbon );
    } else {
      pmtCalib->setPe_per_mip( 0. );
    }
  }

  if ( id.tile() ) {
    pmtCalib->setMip_per_MeV( m_mip_per_MeV );
  } else if ( id.ribbon() ) {
    pmtCalib->setMip_per_MeV( m_mip_per_MeV_ribbon );
  } else {
    pmtCalib->setMip_per_MeV( 0. );
  }

  sc = pmtCalib->latchPePerMeV();

  sc = m_calibSvc->getPedestal(id,pmt,ped);
  if ( sc.isFailure() ) return sc;
  pmtCalib->setPedestal(*ped);

  //Change made by D. Green to incorporate zero suppression as a JO
  sc = pmtCalib->latchPedestal();
  if ( sc.isFailure() ) return sc;
  sc = pmtCalib->latchPhaThreshold(m_counts_above_pedestal);
  if ( sc.isFailure() ) return sc;

  sc = m_calibSvc->getMipPeak(id,pmt,gain);
  if ( sc.isFailure() ) return sc;
  pmtCalib->setMipPeak(*gain);

  sc = m_calibSvc->getHighRange(id,pmt,highRange);
  if ( sc.isFailure() ) return sc;  
  pmtCalib->setHighRange(*highRange);

  sc = m_calibSvc->getRange(id,pmt,range);
  if ( sc.isFailure() ) return sc;  
  pmtCalib->setRangeXOver(*range);
  sc = pmtCalib->latchXOverMips( );
  if ( sc.isFailure() ) return sc;
 
  if ( m_veto_mips > 0. ) {
    pmtCalib->setVetoThresholdMips(m_veto_mips,m_veto_width_mips);
  } else {
    sc = m_calibSvc->getVeto(id,pmt,veto);
    if ( sc.isFailure() ) return sc;
    pmtCalib->setVetoThresh(*veto);
    pmtCalib->setVetoWidthMips(m_veto_width_mips);
    sc = pmtCalib->latchVetoThreshold( );
    if ( sc.isFailure() ) return sc;
  }

  if ( m_cno_mips > 0. ) {
    pmtCalib->setCnoThresholdMips(m_cno_mips,m_cno_width_mips);
  } else {
    sc = m_calibSvc->getCno(id,pmt,cno);
    if ( sc.isFailure() ) return sc;
    pmtCalib->setCnoThresh(*cno);  
    pmtCalib->setCnoWidthMips(m_cno_width_mips);
    sc = pmtCalib->latchCnoThreshold( );
    if ( sc.isFailure() ) return sc;
  }

  sc = m_calibSvc->getCoherentNoise(id,pmt,coherentNoise);
  if ( sc.isFailure() ) return sc;
  pmtCalib->setCoherentNoise(*coherentNoise);    

  if ( id.ribbon() ) {
    sc = m_calibSvc->getRibbon(id,pmt,ribbon);
    if ( sc.isFailure() ) return sc;
    pmtCalib->setRibbon(*ribbon);    
  }

  return sc;
}



