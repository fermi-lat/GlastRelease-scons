#include "../AcdUtil/AcdCalib.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"
#include "CalibData/Acd/AcdRibbon.h"

#include <iostream>

namespace {
    // The following define the standard calib. values
    CalibData::AcdCalibObj::STATUS ok( CalibData::AcdCalibObj::OK );
    CalibData::AcdPed idealPed(       CalibData::AcdPed(0.,0.,ok)        ); // all pedestals are null
    CalibData::AcdGain ribbonGain(    CalibData::AcdGain(56.875,25.4,ok) ); // ribbons
    CalibData::AcdGain tileGain (     CalibData::AcdGain(204.75,50.,ok)  ); // most tiles  
    CalibData::AcdGain tile_12mmGain( CalibData::AcdGain(245.7,50.,ok)   ); // 12mm thick tiles
    CalibData::AcdGain naGain (       CalibData::AcdGain(-1.,0.,ok)      ); // NA channels      
    CalibData::AcdVeto idealVeto (    CalibData::AcdVeto(-1.,0.,ok)      );// veto fires at 50 counts PHA
    CalibData::AcdCno  idealCno (     CalibData::AcdCno(50.,0.,ok)       ); // cno pedestals are ideal 
    CalibData::AcdRange idealRange (  CalibData::AcdRange(4000.,40.,ok));;   // Switch occurs at 4000 in low range = 0 in High Range 
    CalibData::AcdHighRange idealHighRange (CalibData::AcdHighRange(0.,2.04,4000.,ok));// Pedestal = 0, slope = 2.4 PHA/mip, saturates at 4000 PHA
    CalibData::AcdCoherentNoise idealCoherentNoise (CalibData::AcdCoherentNoise(0.,0.,0.,0.,ok));// Amplitude is 0, no oscillation
  
    // for testing (30 PHA counts, time constant of 500 ticks, oscillation of 500 ticks, phase -1.
    CalibData::AcdCoherentNoise realCoherentNoise(CalibData::AcdCoherentNoise(30.,800.,0.0054,-1.,ok)); 
    
    // PMT is on + side of detector
    CalibData::AcdRibbon idealRibbon_Plus( CalibData::AcdRibbon(0.4,0.6,0.8,1.4,2.2,3.0,200.,ok));
    // PMT is on - side of detector
    CalibData::AcdRibbon idealRibbon_Minus( CalibData::AcdRibbon(3.0,2.2,1.4,0.8,0.6,0.4,200.,ok));
}

CalibData::AcdCalibObj* AcdCalib::getIdeal(AcdCalibData::CALTYPE cType, 
                                           idents::AcdId id, unsigned pmt) 
{

    switch ( cType ) {
    case AcdCalibData::PEDESTAL: return &idealPed;
    case AcdCalibData::GAIN: 
    case AcdCalibData::RIBBON: break;
    case AcdCalibData::VETO: return &idealVeto ;
    case AcdCalibData::CNO: return  &idealCno;
    case AcdCalibData::RANGE: return  &idealRange;
    case AcdCalibData::HIGH_RANGE: return &idealHighRange ;
      // switch to test
      //case AcdCalibData::COHERENT_NOISE: return &realCoherentNoise ;     
    case AcdCalibData::COHERENT_NOISE: return &idealCoherentNoise ; 
    default:
        return 0;
    }

    if ( cType == AcdCalibData::GAIN ) {
      if ( id.ribbon() ) {
        return &ribbonGain;
      } else if ( id.tile() ) {
        if ( id.face() == 0 && id.row() == 2 ) {
	  return &tile_12mmGain;
        } else {
	  return &tileGain;
        }
      } 
      return &naGain;
    } else if ( cType == AcdCalibData::RIBBON ) {
      switch ( id.id() ) {
      case 500:
      case 501:
      case 601:
      case 603:
	return pmt == 0 ? &idealRibbon_Plus : &idealRibbon_Minus;
      case 502:
      case 503:
      case 600:
      case 602:
	return pmt == 0 ? &idealRibbon_Minus : &idealRibbon_Plus;
      default:
	return 0;
      }
    }
    return 0;
}


ICalibPathSvc::CalibItem AcdCalib::calibItem(AcdCalibData::CALTYPE cType) 
{
    switch ( cType ) {
    case AcdCalibData::PEDESTAL: return ICalibPathSvc::Calib_ACD_Ped ;
    case AcdCalibData::GAIN: return ICalibPathSvc::Calib_ACD_ElecGain ;
    case AcdCalibData::VETO: return ICalibPathSvc::Calib_ACD_ThreshVeto ;
    case AcdCalibData::CNO: return ICalibPathSvc::Calib_ACD_ThreshHigh ;
    case AcdCalibData::RANGE: return ICalibPathSvc::Calib_ACD_Range ;
    case AcdCalibData::HIGH_RANGE: return ICalibPathSvc::Calib_ACD_HighRange ;
    case AcdCalibData::COHERENT_NOISE: return ICalibPathSvc::Calib_ACD_CoherentNoise ;
    case AcdCalibData::RIBBON: return ICalibPathSvc::Calib_ACD_Ribbon ;
    default:
        ;
    }
    return ICalibPathSvc::Calib_COUNT;
}

