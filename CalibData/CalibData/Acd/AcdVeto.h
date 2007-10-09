// $Header$
#ifndef CalibData_AcdVeto_h
#define CalibData_AcdVeto_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  class AcdVetoFitDesc : public AcdCalibDescription {
  public:    
    static const AcdVetoFitDesc& instance() {
      static const AcdVetoFitDesc desc;
      return desc;
    };     
  public:
    virtual ~AcdVetoFitDesc(){;};    
  private:    
    AcdVetoFitDesc()
      :AcdCalibDescription(AcdCalibData::VETO,"ACD_Veto"){
      addVarName("veto");
      addVarName("width");
    }
  };

  class AcdVeto : public AcdCalibObj {
  public:
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_ThreshVeto;
    }
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::VETO;
    }
  public:
    AcdVeto(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    AcdVeto(float veto, float width, STATUS status) :
      AcdCalibObj(status,AcdVetoFitDesc::instance()){
      setVals(veto,width,status);
    }
    virtual ~AcdVeto() {}

    float getVeto() const { return (*this)[0];}
    float getWidth() const { return (*this)[1]; }
  };
}


#endif
