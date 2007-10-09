// $Header$
#ifndef CalibData_AcdHighRange_h
#define CalibData_AcdHighRange_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  class AcdHighRangeFitDesc : public AcdCalibDescription {
  public:    
    static const AcdHighRangeFitDesc& instance() {
      static const AcdHighRangeFitDesc desc;
      return desc;
    }       
  public:
    virtual ~AcdHighRangeFitDesc(){;};    
  private:    
    AcdHighRangeFitDesc()
      :AcdCalibDescription(AcdCalibData::HIGH_RANGE,"ACD_HighRange"){
      addVarName("pedestal");
      addVarName("slope");
      addVarName("saturation");      
    }
  };

  class AcdHighRange : public AcdCalibObj {
  public:
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_HighRange;
    }
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::HIGH_RANGE;
    }
  public:
    AcdHighRange(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    AcdHighRange(float pedestal, float slope, float saturation, STATUS status) :
      AcdCalibObj(status,AcdHighRangeFitDesc::instance()){
      setVals(pedestal,slope,saturation,status);
    }
    virtual ~AcdHighRange() {}

    float getPedestal() const { return (*this)[0];}
    float getSlope() const { return (*this)[1]; }
    float getSaturation() const { return (*this)[2]; }

  };
}


#endif
