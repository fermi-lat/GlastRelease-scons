
#ifndef __AcdPha2MipTool_H
#define __AcdPha2MipTool_H 1

#include "AcdIPha2MipTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"

#include "AcdUtil/IAcdCalibSvc.h"

class MsgStream;
class HepMatrix;

/**   
* @class AcdPha2MipTool
*
* @brief Gaudi tool which converts raw PHA counts into MIP equivalent signals
*
* The actual conversions live in AcdUtil::AcdCalibFuncs
* This code just loops of there objects calls the relvent conversions routines
* and makes new objects with the new calibrated quantities
*
* TDS inputs:
*  - const Event::AcdDigiCol&, reference passed in call to makeAcdHits()
*    - Collection of all the AcdDigi (ie raw signals) in the event.
*
* TDS outputs:
*  - Event::AcdHitCol&, reference passed in call to makeAcdHits()
*    - Collection of all the AcdHit (ie calibrated signals) in the event.
*
* Algorithm:
*  - For each AcdDigi
*    -# Grab all the flags out of the digi
*    -# Convert to mips in getCalibratedValues()
*       - LOW range, uses getValues_lowRange()
*          - Get the Pedestal and Gain (aka. MIP peak) calibrations for PMT
*          - MIPs = (PHA - pedestal) / mip_peak
*       - HIGH range, uses getValues_highRange()
*          - Get the HighRange calibration for PMT
*          - MIPs = (PHA - pedestal) * saturation * slope / ( saturation + (PHA - pedestal) )
*    -# Check to see if hit passes cuts in accept()
*    -# Build Event::AcdHit object from flags and calibrated values, add it to output collection
*
* This tool has 5 JO paramters:  
*  - AcdCalibSvc ["AcdCalibSvc"]  : Name of Acd Calibration SVC to use
*  - PHATileCut [0.]              : Ignore all tiles with pedestal subtracted PHA below this value
*  - MIPSTileCut [0.]             : Ignore all tiles with MIP equivalent below this value
*  - PHARibbonCut [0.]            : Ignore all ribbons with pedestal subtracted PHA below this valu
*  - MIPSRibbonCut [0.]           : Ignore all ribbons with MIP equivalent below this value
*
* $Header$
*/

class AcdPha2MipTool : public AcdIPha2MipTool,  public AlgTool {
	
public:
  
  /// Standard tool initialization 
  AcdPha2MipTool
  ( const std::string & type, 
    const std::string & name,
    const IInterface * parent ) ;

  /// D'tor is just cleanup
  virtual ~AcdPha2MipTool() ;
  
  /// @brief Intialization of the tool
  virtual StatusCode initialize();
  
  /**
   * @brief Make collection of hits and fill the hit map
   * 
   * @param digiCol all the input AcdDigis
   * @param periodicEvent is the event a periodic trigger or not
   * @param hitCol collection to be filled with output AcdHit objects
   * @param hitMap hap to be filled with output AcdHit objects, maped by AcdId
   * @return Success or Failure
  **/
  virtual StatusCode makeAcdHits ( const Event::AcdDigiCol&,
				   bool periodicEvent,
				   unsigned gemDeltaEventTime, 
				   Event::AcdHitCol&,
				   AcdRecon::AcdHitMap&);
  
  /**
   * @brief Make a single hit
   * 
   * @param digi input AcdDigi
   * @param periodicEvent is the event a periodic trigger or not
   * @param hit newly made output AcdHit object
   * @return Success or Failure
  **/
  virtual StatusCode makeAcdHit ( const Event::AcdDigi&,
				  bool periodicEvent, 
				  Event::AcdHit*& );

protected:
  
  /**
   * @brief 
   * 
   * @param digi input AcdDigi
   * @param mipsPmtA signal in PMT A, in mips
   * @param mipsPmtB signal in PMT B, in mips
   * @param acceptDigi is hit should be accepted, false otherwise
   * @return Success or Failure
  **/  
  bool getCalibratedValues(const Event::AcdDigi& digi, double& mipsPmtA, double& mipsPmtB, bool& acceptDigi) const;

  /**
   * @brief
   * 
   * @param id of channel in question
   * @param pmt A or B side PMT
   * @param pha raw signal 
   * @param pedSub pedestal subtracted signal
   * @param mips signal in mip equivalent
   * @return Success or Failure
  **/
  bool getValues_lowRange(const idents::AcdId& id, Event::AcdDigi::PmtId pmt, unsigned short pha, 
			  double& pedSub, double& mips) const;
  /**
   * @brief
   *
   * @param id of channel in question
   * @param pmt A or B side PMT
   * @param pha raw signal 
   * @param pedSub pedestal subtracted signal
   * @param mips signal in mip equivalent
   * @return Success or Failure
  **/
  bool getValues_highRange(const idents::AcdId& id, Event::AcdDigi::PmtId pmt, unsigned short pha, 
			   double& pedSub, double& mips) const;  

  /**
   * @brief 
   * 
   * @param id of channel in question
   * @param pedSubtracted pedestal subtracted signal
   * @param mips signal in mip equivalent
   * @return true if hit should be accepted, false otherwise
  **/
  bool accept(const idents::AcdId& id, float pedSubtracted, float mips) const;


  /**
   * @brief 
   * 
   * @param id of channel in question
   * @param total mips
   * @return true if hit should be accepted, false otherwise
  **/
  bool acceptTotalMips(const idents::AcdId& id, float totalMips) const;
  

private:

  /// cut tiles with pedestal subtracted PHA below this value
  float m_pha_tile_cut;
  /// cut tiles with MIP equivalent signal below this value
  float m_mips_tile_cut;
  /// cut ribbons with pedestal subtracted PHA below this value
  float m_pha_ribbon_cut;
  /// cut ribbons with MIP equivalent signal below this value
  float m_mips_ribbon_cut;
  /// Value to use for "Ninja" hits, with no signal, but veto bit asserted
  float m_vetoThreshold;
  /// Flag to apply the coherent noise calibration
  bool m_applyCoherentNoiseCalib;

  /// Number of ticks since last readout.  Needed to calibration coherent noise
  unsigned m_gemDeltaEventTime;

  /// Output collection
  Event::AcdHitCol* output;

  /// Name of AcdCalibSvc to get calibrations
  std::string  m_calibSvcName;

  /// AcdCalibSvc to get calibrations
  AcdUtil::IAcdCalibSvc* m_calibSvc;
  

} ;

#endif
	
	
