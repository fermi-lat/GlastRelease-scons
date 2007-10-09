// $Header$
#ifndef CalibData_AcdRange_h
#define CalibData_AcdRange_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  class AcdRangeFitDesc : public AcdCalibDescription {
  public:    
    static const AcdRangeFitDesc& instance() {
      static const AcdRangeFitDesc desc;
      return desc;
    };        
  public:
    virtual ~AcdRangeFitDesc(){;};    
  private:    
    AcdRangeFitDesc()
      :AcdCalibDescription(AcdCalibData::RANGE,"ACD_Range"){
      addVarName("low_max");
      addVarName("high_min");
    }
  };

  class AcdRange : public AcdCalibObj {
  public:
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_Range;
    } 
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::RANGE;
    }
  public:
    AcdRange(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    AcdRange(float lowMax, float highMin, STATUS status) :
      AcdCalibObj(status,AcdRangeFitDesc::instance()){
      setVals(lowMax,highMin,status);
    }
    virtual ~AcdRange() {}

    float getLowMax() const { return (*this)[0];}
    float getHighMin() const { return (*this)[1]; }
  };
}


#endif
