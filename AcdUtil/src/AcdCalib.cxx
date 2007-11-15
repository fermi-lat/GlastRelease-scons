#include "../AcdUtil/AcdCalib.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"

namespace AcdCalib {

  // The following will be used to keep a copy of standard calib. values
  static CalibData::AcdPed* idealPed = 0; 
  static CalibData::AcdGain* ribbonGain = 0; 
  static CalibData::AcdGain* tileGain = 0;   
  static CalibData::AcdGain* tile_12mmGain = 0;
  static CalibData::AcdGain* naGain = 0;       
  static CalibData::AcdVeto* idealVeto = 0;
  static CalibData::AcdCno*  idealCno = 0; 
  static CalibData::AcdRange* idealRange = 0;    
  static CalibData::AcdHighRange* idealHighRange = 0;
  static CalibData::AcdCoherentNoise* idealCoherentNoise = 0;

  CalibData::AcdCalibObj* getIdeal(AcdCalibData::CALTYPE cType, 
                                   idents::AcdId id, unsigned /* pmt */) {

    static CalibData::AcdCalibObj::STATUS ok = CalibData::AcdCalibObj::OK;

    if (idealPed == 0) {
      idealPed = new CalibData::AcdPed(0.,0.,ok); // all pedestals are null

      ribbonGain = new CalibData::AcdGain(56.875,25.4,ok);   // ribbons
      tileGain = new CalibData::AcdGain(204.75,50.,ok);      // most tiles
      tile_12mmGain = new CalibData::AcdGain(245.7,50.,ok); // 12mm thick tiles
      naGain = new CalibData::AcdGain(-1.,0.,ok);            // NA channels

      // veto fires at 50 counts PHA
      idealVeto = new CalibData::AcdVeto(50.,0.,ok); 
      idealCno = new CalibData::AcdCno(50.,0.,ok);  // cno pedestals are ideal

      // Switch occurs at 4000 in low range = 0 in High Range
      idealRange = new CalibData::AcdRange(4000.,40.,ok);
      // Pedestal = 0, slope = 2.4 PHA/mip, saturates at 4000 PHA
      idealHighRange = new CalibData::AcdHighRange(0.,2.04,4000.,ok); 

      // Amplitude is 0, no oscillation
      idealCoherentNoise = new CalibData::AcdCoherentNoise(0.,0.,0.,0.,ok); 
    }
    switch ( cType ) {
    case AcdCalibData::PEDESTAL: return idealPed;
    case AcdCalibData::GAIN: break;
    case AcdCalibData::VETO: return idealVeto ;
    case AcdCalibData::CNO: return  idealCno;
    case AcdCalibData::RANGE: return  idealRange;
    case AcdCalibData::HIGH_RANGE: return idealHighRange ;
    case AcdCalibData::COHERENT_NOISE: return idealCoherentNoise ;  
    default:
      return 0;
    }

    if ( id.ribbon() ) {
      return ribbonGain;
    } else if ( id.tile() ) {
      if ( id.face() == 0 && id.row() == 2 ) {
	return tile_12mmGain;
      } else {
	return tileGain;
      }
    } 
    return naGain;
  }

  ICalibPathSvc::CalibItem calibItem(AcdCalibData::CALTYPE cType) {
    switch ( cType ) {
    case AcdCalibData::PEDESTAL: return ICalibPathSvc::Calib_ACD_Ped ;
    case AcdCalibData::GAIN: return ICalibPathSvc::Calib_ACD_ElecGain ;
    case AcdCalibData::VETO: return ICalibPathSvc::Calib_ACD_ThreshVeto ;
    case AcdCalibData::CNO: return ICalibPathSvc::Calib_ACD_ThreshHigh ;
    case AcdCalibData::RANGE: return ICalibPathSvc::Calib_ACD_Range ;
    case AcdCalibData::HIGH_RANGE: return ICalibPathSvc::Calib_ACD_HighRange ;
    case AcdCalibData::COHERENT_NOISE: return ICalibPathSvc::Calib_ACD_CoherentNoise ;
    default:
      ;
    }
    return ICalibPathSvc::Calib_COUNT;
  }

}
