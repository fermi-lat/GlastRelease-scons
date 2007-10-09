// $Header$
#ifndef CalibData_AcdGain_h
#define CalibData_AcdGain_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  class AcdGainFitDesc : public AcdCalibDescription {
  public:
    static const AcdGainFitDesc & instance() {
      static const AcdGainFitDesc desc;
      return desc;
    }
  public:
    virtual ~AcdGainFitDesc(){;};    
  private:    
    AcdGainFitDesc()
      :AcdCalibDescription(AcdCalibData::GAIN,"ACD_Gain"){
      addVarName("peak");
      addVarName("width");
    }
  };

  class AcdGain : public AcdCalibObj {    
  public:
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_ElecGain;
    }
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::GAIN;
    }
  public:
    AcdGain(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    AcdGain(float peak, float width, STATUS status) :
      AcdCalibObj(status,AcdGainFitDesc::instance()){
      setVals(peak,width,status);
    }
    virtual ~AcdGain() {}

    float getPeak() const { return (*this)[0];}
    float getWidth() const { return (*this)[1]; }
  };
}

#endif
