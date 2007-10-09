#ifndef CalibData_AcdCalib_h
#define CalibData_AcdCalib_h

#include "CalibData/Acd/AcdCalibBase.h"
#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"

namespace CalibData {

  template <class CalibObjType> 
  class AcdCalib : public AcdCalibBase {

  public:

    // this has to be here so that other template know what 
    // type of object is managed by this class
    typedef CalibObjType ObjType;

    static AcdCalibData::CALTYPE calibType()  {
      return CalibObjType::calibType();
    }
    
  public:

    AcdCalib(const AcdCalibDescription& desc,
	     unsigned nFace=5, unsigned nRow=5, unsigned nCol=5, 
	     unsigned nNA=11, unsigned nPmt=2)
      :AcdCalibBase(desc,nFace,nRow,nCol,nNA,nPmt){;}

    virtual ~AcdCalib() {;}

    virtual const CLID& clID() const {
      return CalibObjType::calibCLID();
    }

    
    CalibObjType* getPmt(idents::AcdId id,unsigned pmt) {
      return static_cast<CalibObjType*>(get(id,pmt));
    }

    bool putPmt(idents::AcdId id,unsigned pmt, CalibObjType& pmtCalib) {
      return put(id,pmt,pmtCalib);
    }

  protected:

    virtual AcdCalibObj* makeNew() const {
      static std::vector<float> nullVect;
      return new CalibObjType(*(desc()),nullVect);
    }
  };

  typedef AcdCalib<CalibData::AcdPed> AcdPedCalib;  
  typedef AcdCalib<CalibData::AcdGain> AcdGainCalib;
  typedef AcdCalib<CalibData::AcdVeto> AcdVetoCalib;
  typedef AcdCalib<CalibData::AcdCno> AcdCnoCalib;
  typedef AcdCalib<CalibData::AcdRange> AcdRangeCalib; 
  typedef AcdCalib<CalibData::AcdHighRange> AcdHighRangeCalib;  
  typedef AcdCalib<CalibData::AcdCoherentNoise> AcdCoherentNoiseCalib;

}

#endif
