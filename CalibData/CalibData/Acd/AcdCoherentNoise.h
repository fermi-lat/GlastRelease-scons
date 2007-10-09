// $Header$
#ifndef CalibData_AcdCoherentNoise_h
#define CalibData_AcdCoherentNoise_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  class AcdCoherentNoiseFitDesc : public AcdCalibDescription {
  public:    
    static const AcdCoherentNoiseFitDesc& instance(){
      static const AcdCoherentNoiseFitDesc desc;
      return desc;
    }      
  public:
    virtual ~AcdCoherentNoiseFitDesc(){;};    
  private:    
    AcdCoherentNoiseFitDesc()
      :AcdCalibDescription(AcdCalibData::COHERENT_NOISE,"ACD_CoherentNoise"){
      addVarName("veto");
      addVarName("width");
    }
  };


  class AcdCoherentNoise : public AcdCalibObj {
  public:
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_CoherentNoise;
    }
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::COHERENT_NOISE;
    }
  public:
    AcdCoherentNoise(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    AcdCoherentNoise(float amplitude, float decay, float frequency, float phase, STATUS status) :
      AcdCalibObj(status,AcdCoherentNoiseFitDesc::instance()){
      setVals(amplitude,decay,frequency,phase,status);
    }
    virtual ~AcdCoherentNoise() {}
    float getAmplitude() const { return (*this)[0];}
    float getDecay() const { return (*this)[1]; }
    float getFrequency() const { return (*this)[2]; }
    float getPhase() const { return (*this)[3]; }
  };
}

#endif
