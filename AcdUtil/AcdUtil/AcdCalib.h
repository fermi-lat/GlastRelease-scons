#ifndef ACDCALIB_H
#define ACDCALIB_H

#include "CalibData/Acd/AcdCalibEnum.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "idents/AcdId.h" 

namespace CalibData {
  class AcdCalibObj;
};

namespace AcdCalib {
  


  CalibData::AcdCalibObj* getIdeal(AcdCalibData::CALTYPE cType, idents::AcdId id, unsigned pmt);

  ICalibPathSvc::CalibItem calibItem(AcdCalibData::CALTYPE cType); 

}

#endif
