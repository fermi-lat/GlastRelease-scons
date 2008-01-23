#ifndef ACDCALIB_H
#define ACDCALIB_H

#include "CalibData/Acd/AcdCalibEnum.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "idents/AcdId.h" 

// Forward declarations
namespace CalibData {
  class AcdCalibObj;
};


/**
 * @brief Some functions used in handling ACD calibrations
 *
 */

namespace AcdCalib {
  

  /**
   * @brief Get the ideal calibration of a given type for one PMT
   *
   * @param cType  an enum with the calibration type
   * @param id     the ID of the pmt in question
   * @param pmt    the PMT index (0 for A side or 1 for B side) 
   * @return a pointer to an object that encapsultes the "ideal" values of the calibration
   */
  CalibData::AcdCalibObj* getIdeal(AcdCalibData::CALTYPE cType, idents::AcdId id, unsigned pmt);

  /**
   * @brief Get the set of all calibrations of a given type
   *
   * @param cType  an enum with the calibration type
   * @return a pointer to an object that encapsultes all the calibrations
   */  
  ICalibPathSvc::CalibItem calibItem(AcdCalibData::CALTYPE cType); 

}

#endif
