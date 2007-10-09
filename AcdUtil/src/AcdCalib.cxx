#include "../AcdUtil/AcdCalib.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"

namespace AcdCalib {


  CalibData::AcdCalibObj* getIdeal(AcdCalibData::CALTYPE cType, idents::AcdId id, unsigned /* pmt */) {

    static CalibData::AcdCalibObj::STATUS ok = CalibData::AcdCalibObj::OK;

    static CalibData::AcdPed idealPed(0.,0.,ok); // all pedestals are null

    // four kinds of channels for Gain
    static CalibData::AcdGain ribbonGain(56.875,25.4,ok);   // ribbons
    static CalibData::AcdGain tileGain(204.75,50.,ok);      // most tiles
    static CalibData::AcdGain tile_12mmGain(245.7,50.,ok);  // 12mm thick tiles
    static CalibData::AcdGain naGain(-1.,0.,ok);            // NA channels

    static CalibData::AcdVeto idealVeto(50.,0.,ok);          // veto fires at 50 counts PHA
    static CalibData::AcdCno  idealCno(50.,0.,ok);           // cno pedestals are ideal

    static CalibData::AcdRange idealRange(4000.,40.,ok);             // Switch occurs at 4000 in low range = 0 in High Range
    static CalibData::AcdHighRange idealHighRange(0.,2.04,4000.,ok); // Pedestal = 0, slope = 2.4 PHA/mip, saturates at 4000 PHA

    static CalibData::AcdCoherentNoise idealCoherentNoise(0.,0.,0.,0.,ok); // Amplitude is 0, no oscillation

    switch ( cType ) {
    case AcdCalibData::PEDESTAL: return &idealPed;
    case AcdCalibData::GAIN: break;
    case AcdCalibData::VETO: return &idealVeto ;
    case AcdCalibData::CNO: return  &idealCno;
    case AcdCalibData::RANGE: return  &idealRange;
    case AcdCalibData::HIGH_RANGE: return &idealHighRange ;
    case AcdCalibData::COHERENT_NOISE: return &idealCoherentNoise ;  
    default:
      return 0;
    }

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
