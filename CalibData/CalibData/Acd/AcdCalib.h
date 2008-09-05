#ifndef CalibData_AcdCalib_h
#define CalibData_AcdCalib_h

/** 
 * @class CalibData::AcdCalib<CalibObjType>
 *
 * @brief Calibrations of a particular type for all the ACD channels.
 *
 * This template provides type safety for the various calibrations.
 * While CalibData::AcdCalibBase gives access to individual channel calibrations, t
 * this class recasts those objects to the correct type.
 * 
 * @author Eric Charles
 * $Header$
 */


#include "CalibData/Acd/AcdCalibBase.h"
#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"
#include "CalibData/Acd/AcdRibbon.h"
#include "CalibData/Acd/AcdPE.h"

namespace CalibData {

  template <class CalibObjType> 
  class AcdCalib : public AcdCalibBase {

  public:

    /// this has to be here so that other template know what 
    /// type of object is managed by this class
    typedef CalibObjType ObjType;

    /// Return the type of calibration
    static AcdCalibData::CALTYPE calibType()  {
      return CalibObjType::calibType();
    }
    
  public:

    /// Allocate space for the calibration curves
    AcdCalib(const AcdCalibDescription& desc,
	     unsigned nFace=5, unsigned nRow=5, unsigned nCol=5, 
	     unsigned nNA=11, unsigned nPmt=2)
      :AcdCalibBase(desc,nFace,nRow,nCol,nNA,nPmt){;}
    
    /// Trivial d'tor
    virtual ~AcdCalib() {;}

    /// For Gaudi
    virtual const CLID& clID() const {
      return CalibObjType::calibCLID();
    }

    /// Get a calibration, cast to the correct type
    CalibObjType* getPmt(idents::AcdId id,unsigned pmt) {
      return static_cast<CalibObjType*>(get(id,pmt));
    }

    /// Add a calibration, 
    bool putPmt(idents::AcdId id,unsigned pmt, CalibObjType& pmtCalib) {
      return put(id,pmt,pmtCalib);
    }

  protected:

    /// Allocate a new calibration for a single channel
    virtual AcdCalibObj* makeNew() const {
      static std::vector<float> nullVect;
      return new CalibObjType(*(desc()),nullVect);
    }
  };

  // declare all the calibration types
  typedef AcdCalib<CalibData::AcdPed> AcdPedCalib;  
  typedef AcdCalib<CalibData::AcdGain> AcdGainCalib;
  typedef AcdCalib<CalibData::AcdVeto> AcdVetoCalib;
  typedef AcdCalib<CalibData::AcdCno> AcdCnoCalib;
  typedef AcdCalib<CalibData::AcdRange> AcdRangeCalib; 
  typedef AcdCalib<CalibData::AcdHighRange> AcdHighRangeCalib;  
  typedef AcdCalib<CalibData::AcdCoherentNoise> AcdCoherentNoiseCalib;
  typedef AcdCalib<CalibData::AcdRibbon> AcdRibbonCalib;
  typedef AcdCalib<CalibData::AcdPE> AcdPECalib;

}

#endif
