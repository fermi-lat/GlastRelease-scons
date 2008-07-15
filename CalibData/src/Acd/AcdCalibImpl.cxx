// $Header$


#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"
#include "CalibData/Acd/AcdRibbon.h"
#include "CalibData/Acd/AcdHighPed.h"
#include "CalibData/Acd/AcdCarbon.h"
#include "CalibData/Acd/AcdVetoFit.h"
#include "CalibData/Acd/AcdCnoFit.h"

namespace CalibData {

  int buildAcdCalibDescriptions() {
    // careful, order matters within a calibration type

    // pedestals
    AcdPedestalFitDesc::instance();
    // gains (aka mip peaks)
    AcdGainFitDesc::instance();
    // veto set points 
    AcdVetoFitDesc::instance();
    // cno set points
    AcdCnoFitDesc::instance();
    // range crossover
    AcdRangeFitDesc::instance();
    // high range calibration
    AcdHighRangeFitDesc::instance();
    // coherent noise 
    AcdCoherentNoiseFitDesc::instance();    
    // ribbons
    AcdRibbonFitDesc::instance();    
    // high range pedestals
    AcdHighPedestalFitDesc::instance();   
    // Carbon peak calibrations
    AcdCarbonFitDesc::instance();
    // Veto Fit parameters
    AcdVetoFitFitDesc::instance();
    // CNO Fit parameters
    AcdCnoFitFitDesc::instance();
    
    return 0;
  }
    
}

int dummy = CalibData::buildAcdCalibDescriptions();
