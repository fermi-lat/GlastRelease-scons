// $Header$
#ifndef CalibData_AcdCoherentNoise_h
#define CalibData_AcdCoherentNoise_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdCoherentNoiseFitDesc
   *
   * @brief Description of coherent readout noise calibration.
   * 
   * This calibration consists of:
   *  - amplitude   = the amplitude of the coherent noise (in PHA counts)
   *  - decay       = the decay constant of the coherent noise (in 50MHz ticks)
   *  - frequency   = the frequency of the coherent noise (in 1./50MHz ticks)
   *  - phase       = the phase of the coherent noise (in 50MHz ticks)
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCoherentNoiseFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdCoherentNoiseFitDesc& instance(){
      static const AcdCoherentNoiseFitDesc desc;
      return desc;
    }      
  public:
    /// Trivial D'ctor
    virtual ~AcdCoherentNoiseFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdCoherentNoiseFitDesc()
      :AcdCalibDescription(AcdCalibData::COHERENT_NOISE,"ACD_CoherentNoise"){
      addVarName("amplitude");
      addVarName("decay");
      addVarName("frequency");
      addVarName("phase");
    }
  };

  /** 
   * @class AcdCoherentNoise
   *
   * @brief Coherent readout noise calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - amplitude   = the amplitude of the coherent noise (in PHA counts)
   *  - decay       = the decay constant of the coherent noise (in 50MHz ticks)
   *  - frequency   = the frequency of the coherent noise (in 1./50MHz ticks)
   *  - phase       = the phase of the coherent noise (in 50MHz ticks)
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCoherentNoise : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_CoherentNoise;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::COHERENT_NOISE;
    }
  public:
    /// Build from description and a set of values
    AcdCoherentNoise(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdCoherentNoise(float amplitude, float decay, float frequency, float phase, STATUS status) :
      AcdCalibObj(status,AcdCoherentNoiseFitDesc::instance()){
      setVals(amplitude,decay,frequency,phase,status);
    }
    /// Trivial d'tor
    virtual ~AcdCoherentNoise() {}

    // Provide access to the values
    float getAmplitude() const { return (*this)[0];}
    float getDecay() const { return (*this)[1]; }
    float getFrequency() const { return (*this)[2]; }
    float getPhase() const { return (*this)[3]; }
  };
}

#endif
