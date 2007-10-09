// $Header$


#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"

namespace CalibData {

  int buildAcdCalibDescriptions() {
    // careful, order matters within a calibration type

    // pedestals
    const AcdPedestalFitDesc& ped = AcdPedestalFitDesc::instance();
    // gains (aka mip peaks)
    const AcdGainFitDesc& gain = AcdGainFitDesc::instance();
    // veto set points 
    const AcdVetoFitDesc& veto = AcdVetoFitDesc::instance();
    // cno set points
    const AcdCnoFitDesc& cno = AcdCnoFitDesc::instance();
    // range crossover
    const AcdRangeFitDesc& range = AcdRangeFitDesc::instance();
    // high range calibration
    const AcdHighRangeFitDesc& highRange = AcdHighRangeFitDesc::instance();
    // coherent noise 
    const AcdCoherentNoiseFitDesc& coherentNois = AcdCoherentNoiseFitDesc::instance();    
    return 0;
  }
    
}

int dummy = CalibData::buildAcdCalibDescriptions();
