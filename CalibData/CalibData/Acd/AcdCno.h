// $Header$
#ifndef CalibData_AcdCno_h
#define CalibData_AcdCno_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  class AcdCnoFitDesc : public AcdCalibDescription {
  public:    
    static const AcdCnoFitDesc& instance() {
      static const AcdCnoFitDesc desc;
      return desc;
    };        
  public:
    virtual ~AcdCnoFitDesc(){;};    
  private:    
    AcdCnoFitDesc()
      :AcdCalibDescription(AcdCalibData::CNO,"ACD_Cno"){
      addVarName("cno");
      addVarName("width");
    }
  };

  class AcdCno : public AcdCalibObj {
  public:
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_ThreshHigh;
    }
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::CNO;
    }
  public:
    AcdCno(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    AcdCno(float cno, float width, STATUS status) :
      AcdCalibObj(status,AcdCnoFitDesc::instance()){    
      setVals(cno,width,status);
    }
    virtual ~AcdCno() {}

    float getCno() const { return (*this)[0];}
    float getWidth() const { return (*this)[1]; }
  };
}


#endif
