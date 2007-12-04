#include "../AcdUtil/AcdCalib.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"

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
    CalibData::AcdHighRange idealHighRange 
        (CalibData::AcdHighRange(0.,2.04,4000.,ok));// Pedestal = 0, slope = 2.4 PHA/mip, saturates at 4000 PHA
    CalibData::AcdCoherentNoise 
        idealCoherentNoise (CalibData::AcdCoherentNoise(0.,0.,0.,0.,ok));// Amplitude is 0, no oscillation
}

CalibData::AcdCalibObj* AcdCalib::getIdeal(AcdCalibData::CALTYPE cType, 
                                           idents::AcdId id, unsigned /* pmt */) 
{

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
    default:
        ;
    }
    return ICalibPathSvc::Calib_COUNT;
}

