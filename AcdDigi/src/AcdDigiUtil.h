#ifndef _AcdDigi_AcdDigiUtil_H
#define _AcdDigi_AcdDigiUtil_H 1


// stl
#include <map>
#include <string>

// GLAST
#include "GaudiKernel/StatusCode.h"
#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "Event/Digi/AcdDigi.h"


class IAcdGeometrySvc;
class MsgStream;

// Forward Declartions
namespace AcdUtil {
  class IAcdCalibSvc;
};

namespace xmlBase {
  class IFile;
};

namespace CalibData {
  class AcdPed;
  class AcdGain;
  class AcdHighRange;
  class AcdRange;
  class AcdVeto;
  class AcdCno;
  class AcdCoherentNoise;
  class AcdRibbon;
};

namespace Event {
  class McPositionHit;
};

/** 
 * @class AcdSimCalibData
 *
 * @brief Utility class that bundles all calibrations constants associated with a single ACD channel.
 * 
 * @author Eric Charles
 * $Header$
 */

class AcdSimCalibData {
public: 

  /// Trivial c'tor
  AcdSimCalibData();
    
  /// Trivial d'tor
  virtual ~AcdSimCalibData(){;}

  /// photo-electrons per mip
  inline double pe_per_mip() const { return m_pe_per_mip; } 
  /// mip-equivalent light yield by MeV
  inline double mip_per_MeV() const { return m_mip_per_MeV; }
  /// photo-electrons by MeV
  inline double pe_per_MeV() const { return m_pe_per_MeV; }
  /// pedestal
  inline double pedestal_pha() const{return m_pedestal;}
  /// zero suppression threshold
  inline unsigned short threshold_pha() const { return m_threshold_pha; }
  /// veto threshold  
  inline double veto_threshold_mips() const { return m_veto_threshold_mips; }
  inline double veto_width_mips() const { return m_veto_width_mips; }
  /// range crossover
  inline double xover_mips() const { return m_xover_mips; }
  /// cno threshold
  inline double cno_threshold_mips() const { return m_cno_threshold_mips; }
  inline double cno_width_mips() const { return m_cno_width_mips; }
    
  // Lots of access functions for setting calibration objects
  inline void setPedestal(CalibData::AcdPed& ped) { m_ped = &ped; };
  inline void setMipPeak(CalibData::AcdGain& gain) { m_gain = &gain; };
  inline void setHighRange(CalibData::AcdHighRange& highRange) { m_highRange = &highRange; };
  inline void setRangeXOver(CalibData::AcdRange& range) { m_range = &range; };  
  inline void setVetoThresh(CalibData::AcdVeto& veto) { m_veto = &veto; };
  inline void setCnoThresh(CalibData::AcdCno& cno) { m_cno = &cno; };
  inline void setCoherentNoise(CalibData::AcdCoherentNoise& coherentNoise) { m_coherentNoise = &coherentNoise; };
  inline void setRibbon(CalibData::AcdRibbon& ribbon) { m_ribbon = &ribbon; };
  
  // Lots of access functions for setting other quantities
  inline void setPe_per_mip(double val) {  m_pe_per_mip = val; };
  inline void setMip_per_MeV(double val) {  m_mip_per_MeV = val; };
  inline void setPhaThreshold(unsigned short val) {  m_threshold_pha = val; };
  inline void setVetoThresholdMips(double veto, double width) {  m_veto_threshold_mips = veto;  m_veto_width_mips = width; };
  inline void setVetoWidthMips(double width) { m_veto_width_mips = width; };
  inline void setCnoThresholdMips(double cno, double width) {  m_cno_threshold_mips = cno;  m_cno_width_mips = width; };
  inline void setCnoWidthMips(double width) { m_cno_width_mips = width; };
  inline void setXOverMips(double val) {  m_xover_mips = val; };

  // Access to calibrations
  inline CalibData::AcdPed* ped_calib() const { return m_ped; }
  inline CalibData::AcdGain* gain_calib() const { return m_gain; }
  inline CalibData::AcdHighRange* highRange_calib() const { return m_highRange; }
  inline CalibData::AcdRange* range_calib() const { return m_range; }
  inline CalibData::AcdVeto* veto_calib() const { return m_veto; }
  inline CalibData::AcdCno* cno_calib() const { return m_cno; }
  inline CalibData::AcdCoherentNoise* coherentNoise_calib() const { return m_coherentNoise; }
  inline CalibData::AcdRibbon* ribbon_calib() const { return m_ribbon; }

  // Grap combined constants as needed
  StatusCode latchPePerMeV();
  StatusCode latchPhaThreshold(double countsAbovePed);
  //Change made by D. Green to incorporate zero suppression as a JO
  StatusCode latchPedestal();
  StatusCode latchVetoThreshold();  
  StatusCode latchCnoThreshold();  
  StatusCode latchXOverMips();    

private:  

  // calibrations
  CalibData::AcdPed* m_ped;
  CalibData::AcdGain* m_gain;
  CalibData::AcdHighRange* m_highRange;
  CalibData::AcdRange* m_range;  
  CalibData::AcdVeto* m_veto;
  CalibData::AcdCno* m_cno;
  CalibData::AcdCoherentNoise* m_coherentNoise;  
  CalibData::AcdRibbon* m_ribbon;  

  // Derived/ set values

  /// photo-electrons per mip
  double m_pe_per_mip;  
  /// mip-equivalent light yield by MeV
  double m_mip_per_MeV;  
  /// photo-electrons per MeV
  double m_pe_per_MeV;
  /// Zero suppresion threshold
  unsigned short m_threshold_pha;       
  /// veto threshold  
  double m_veto_threshold_mips;
  double m_veto_width_mips;
  /// range crossover
  double m_xover_mips;
  /// cno threshold
  double m_cno_threshold_mips;
  double m_cno_width_mips;
  /// high range calibraion
  double m_pedestal_highRange;
  /// low range pedestal
  double m_pedestal;

};

/** @class AcdDigiUtil
* @brief Utility class that defines the methods used for ACD digitization.
* 
* @author Heather Kelly
* $Header$
*/

class AcdDigiUtil  {

public:

  /// Returns the number of pe seen in a tube
  static double simulateDynodeChain(double pe_meanValue);
  
  /// Returns a value sampled from a Poisson distribution
  /// @param pmtPhotoElectrons is the mean of the Poisson distribution
  static double shootPoisson(double pmtPhotoElectrons);

  /// Returns a value sampled from a Gaussian distribution
  /// @param std_dev Standard Deviation to be  used when sampling
  static double shootGaussian(double std_dev);
  
  /// Compares two volume IDs, returns true if the screw belongs with the tile
  static bool compareVolIds(const idents::VolumeIdentifier& tileId,
			    const idents::VolumeIdentifier& screwId);
  

  /**
   * @brief Return the localY of a ribbon from the MCPositionHit
   * @param hit the MCPositionHit
   * @param geomSvc the geometery service
   * @param ribbonLength length along the ribbon (measured from the center)   
   * @param ribbonBin which bin to use for attenuation
   * @return true for success, false otherwise
   **/
  static bool getRibbonLengthAndBin(const Event::McPositionHit* hit,
				    IAcdGeometrySvc& geomSvc,
				    MsgStream& log,
				    double& ribbonLength, int& ribbonBin);
			  

public:
    
  AcdDigiUtil(); 
  
  virtual ~AcdDigiUtil();

  StatusCode initialize(AcdUtil::IAcdCalibSvc& calibSvc, IAcdGeometrySvc& geomSvc,
			const std::string& xmlFileName);

  /// calulates the number of pe seen in each tube given the 
  /// deposited energy  
  StatusCode photoElectronsFromEnergy(const idents::AcdId& id, double energy,
				      double& pe_pmtA, double& pe_pmtB);

  /// calulates the number of pe seen in each tube given the 
  /// deposited energy
  StatusCode photoElectronsFromEnergy_tile(const Event::McPositionHit *hit, bool edgeEffect, MsgStream& log,
					   double& pe_pmtA, double& pe_pmtB);
  
  /// calulates the number of pe seen in each tube given the 
  /// deposited energy
  StatusCode photoElectronsFromEnergy_ribbon(const Event::McPositionHit *hit, MsgStream& log,
					     double& pe_pmtA, double& pe_pmtB);   
  
  /// calulates the light yield expressed in mip equivalent from the number of observed PE 
  StatusCode mipEquivalentLightYeild(const idents::AcdId& id, double pe_pmtA, double pe_pmtB, MsgStream& log,
				     double& mipEquivA, double& mipEquivB);

  /// get the PHA counts from the mip equivalent light yield
  StatusCode phaCounts(const idents::AcdId& id, const double mipEquiv[2], bool applyNoise, MsgStream& log,
		       Event::AcdDigi::Range range[2], unsigned short pha[2]);

  /// get the PHA counts from the mip equivalent light yield
  StatusCode applyCoherentNoiseToPha(const idents::AcdId& id, unsigned int deltaGemEventTime,  MsgStream& log,
				     unsigned short pha[2]);    

  /// Adjusts the deposited energy recorded in an ACD volume 
  /// based on the location of the hit
  StatusCode tileEdgeEffect(const Event::McPositionHit *hit,  MsgStream& log, double& energy);

  /// Checks all the various thresholds
  //Change made by D. Green to incorporate zero suppression as a JO
  StatusCode checkThresholds(const idents::AcdId& id, const double mipEquiv[2],
			     const unsigned short phaArr[2], const Event::AcdDigi::Range rangeArr[2], bool applyNoise, 
			     MsgStream& log,
			     bool& makeDigi, 
			     bool phaThreshArr[2], bool vetoArr[2], bool highArr[2],
			     const double phaZeroThreshold);  
  
protected:

  // get the calibration data from the local map, of fetch in from the DB if needed
  StatusCode getCalibData(const idents::AcdId& id, AcdSimCalibData*& pmtACalib, AcdSimCalibData*& pmtBCalib);
  
  // fetch the calibration data from the DB
  StatusCode fetchCalibData(const idents::AcdId& id, Event::AcdDigi::PmtId pmt, AcdSimCalibData*& pmtCalib);  

  /// Read data from the input XML file
  StatusCode getParameters(const std::string &xmlFile);

private:
  
  /// The calibration service
  AcdUtil::IAcdCalibSvc* m_calibSvc;
  
  /// The geometry service
  IAcdGeometrySvc* m_geomSvc;

  /// The xml input file
  xmlBase::IFile *m_ifile;
  
  /// Distance (mm) cutoff for applying edge effects 
  double m_max_edge_dist;
  /// Slope of the linear function used to estimate the edge effect
  double m_edge_slope;
  /// y-intercept of the linear function used to estimate the edge effect
  double m_edge_intercept;
  
  /// Global ratio of photoelectrons to mips
  double m_mean_pe_per_mip;
  double m_mean_pe_per_mip_ribbon;  

  /// MIP per MeV
  double m_mip_per_MeV;
  double m_mip_per_MeV_ribbon;

  double m_veto_mips;
  double m_veto_width_mips;
  double m_cno_mips;
  double m_cno_width_mips;

  /// counts above pedestal to set PHA threshold
  double m_counts_above_pedestal;

  /// store AcdId specific number of photoelectrons per mip as they are read in
  std::map< unsigned int, std::pair<AcdSimCalibData*,AcdSimCalibData*> > m_calibMap;

  
};

#endif
