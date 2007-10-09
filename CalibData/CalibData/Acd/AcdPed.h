// $Header$
#ifndef CalibData_AcdPed_h
#define CalibData_AcdPed_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  class AcdPedestalFitDesc : public AcdCalibDescription {
  public:
    static const AcdPedestalFitDesc& instance() {
      static const AcdPedestalFitDesc desc;
      return desc;
    }
  public:
    virtual ~AcdPedestalFitDesc(){;};    
  private:    
    AcdPedestalFitDesc()
      :AcdCalibDescription(AcdCalibData::PEDESTAL,"ACD_Ped"){
      addVarName("mean");
      addVarName("width");
    }
  };


  class AcdPed : public AcdCalibObj {
  public:
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_Ped;
    }
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::PEDESTAL;
    }
  public:
    // Build from a particular descripton
    AcdPed(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    // Description is known
    AcdPed(float mean, float width, STATUS status) :
      AcdCalibObj(status,AcdPedestalFitDesc::instance()){
      setVals(mean,width,status);
    }
    virtual ~AcdPed() {;}
    float getMean() const {return (*this)[0];}
    float getWidth() const {return (*this)[1]; }

  };
}

#endif
