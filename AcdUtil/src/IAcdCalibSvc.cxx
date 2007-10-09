// $Header $
/** @file
    @author Eric Charles
 */
// @file
//
//
// Author: Eric Charles

// LOCAL

#include "AcdUtil/IAcdCalibSvc.h"

#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdCalib.h"


namespace AcdUtil {

  /** \brief get pedestals for given channel
      \param id  the tile or ribbon id
      \param pmt A(0) or B(1) pmt
      \param pedestal a pointer to the relevent pedestal data
 */  
  StatusCode IAcdCalibSvc::getPedestal(idents::AcdId id, unsigned pmt,
				       CalibData::AcdPed*& pedestal) {
    AcdCalibMgr* calibMgr(0);
    StatusCode sc = getCalibMgr(AcdCalibData::PEDESTAL,calibMgr);
    if ( sc.isFailure() ) {
      pedestal = 0;
      return sc;
    }
    AcdCalibMgrTmpl< CalibData::AcdPedCalib >* typeMgr = static_cast< AcdCalibMgrTmpl< CalibData::AcdPedCalib >* >(calibMgr);
    return typeMgr->getCalibration(id,pmt,pedestal);
  };
  
  /** \brief get mip peak for a given channel
      \param id  the tile or ribbon id
      \param pmt A(0) or B(1) pmt
      \param pedestal a pointer to the relevent gain data
  */
  StatusCode IAcdCalibSvc::getMipPeak(idents::AcdId id, unsigned pmt,
				      CalibData::AcdGain*& mipPeak) {
    AcdCalibMgr* calibMgr(0);
    StatusCode sc = getCalibMgr(AcdCalibData::GAIN,calibMgr);
    if ( sc.isFailure() ) {
      mipPeak = 0;
      return sc;
    }
    AcdCalibMgrTmpl< CalibData::AcdGainCalib >* typeMgr = static_cast< AcdCalibMgrTmpl< CalibData::AcdGainCalib >* >(calibMgr);
    return typeMgr->getCalibration(id,pmt,mipPeak);
  };

  /** \brief get veto threshold for a given channel
      \param id  the tile or ribbon id
      \param pmt A(0) or B(1) pmt
      \param veto a pointer to the relevent threshold data
  */
  StatusCode IAcdCalibSvc::getVeto(idents::AcdId id, unsigned pmt,
				   CalibData::AcdVeto*& veto){
    AcdCalibMgr* calibMgr(0);
    StatusCode sc = getCalibMgr(AcdCalibData::VETO,calibMgr);
    if ( sc.isFailure() ) {
      veto = 0;
      return sc;
    }
    AcdCalibMgrTmpl< CalibData::AcdVetoCalib >* typeMgr = static_cast< AcdCalibMgrTmpl< CalibData::AcdVetoCalib >* >(calibMgr);
    return typeMgr->getCalibration(id,pmt,veto);
  };
  
  /** \brief get cno threshold for a given channel
      \param id  the tile or ribbon id
      \param pmt A(0) or B(1) pmt
      \param cno a pointer to the relevent threshold data
  */
  StatusCode IAcdCalibSvc::getCno(idents::AcdId id, unsigned pmt,
				  CalibData::AcdCno*& cno){
    AcdCalibMgr* calibMgr(0);
    StatusCode sc = getCalibMgr(AcdCalibData::CNO,calibMgr);
    if ( sc.isFailure() ) {
      cno = 0;
      return sc;
    }
    AcdCalibMgrTmpl< CalibData::AcdCnoCalib >* typeMgr = static_cast< AcdCalibMgrTmpl< CalibData::AcdCnoCalib >* >(calibMgr);
    return typeMgr->getCalibration(id,pmt,cno);
  };
  
  /** \brief get range crossover for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param range a pointer to the relevent data
    */
  StatusCode IAcdCalibSvc::getRange(idents::AcdId id, unsigned pmt,
				    CalibData::AcdRange*& range){
    AcdCalibMgr* calibMgr(0);
    StatusCode sc = getCalibMgr(AcdCalibData::RANGE,calibMgr);
    if ( sc.isFailure() ) {
      range = 0;
      return sc;
    }
    AcdCalibMgrTmpl< CalibData::AcdRangeCalib >* typeMgr = static_cast< AcdCalibMgrTmpl< CalibData::AcdRangeCalib >* >(calibMgr);
    return typeMgr->getCalibration(id,pmt,range);
  };
  
  /** \brief get high range calibration for a given channel
      \param id  the tile or ribbon id
      \param pmt A(0) or B(1) pmt
      \param highRange a pointer to the relevent high range data
  */
  StatusCode IAcdCalibSvc::getHighRange(idents::AcdId id, unsigned pmt,
					CalibData::AcdHighRange*& highRange){
    AcdCalibMgr* calibMgr(0);
    StatusCode sc = getCalibMgr(AcdCalibData::HIGH_RANGE,calibMgr);
    if ( sc.isFailure() ) {
      highRange = 0;
      return sc;
    }
    AcdCalibMgrTmpl< CalibData::AcdHighRangeCalib >* typeMgr = static_cast< AcdCalibMgrTmpl< CalibData::AcdHighRangeCalib >* >(calibMgr);
    return typeMgr->getCalibration(id,pmt,highRange);
  };
  
  /** \brief get coherent noise calibration for a given channel
      \param id  the tile or ribbon id
      \param pmt A(0) or B(1) pmt
      \param calib a pointer to the relevent high range data
  */
  StatusCode IAcdCalibSvc::getCoherentNoise(idents::AcdId id, unsigned pmt,
					    CalibData::AcdCoherentNoise*& noiseCalib){
    AcdCalibMgr* calibMgr(0);
    StatusCode sc = getCalibMgr(AcdCalibData::COHERENT_NOISE,calibMgr);
    if ( sc.isFailure() ) {
       noiseCalib= 0;
      return sc;
    }
    AcdCalibMgrTmpl< CalibData::AcdCoherentNoiseCalib >* typeMgr = static_cast< AcdCalibMgrTmpl< CalibData::AcdCoherentNoiseCalib >* >(calibMgr);
    return typeMgr->getCalibration(id,pmt,noiseCalib);
  };
  

}
