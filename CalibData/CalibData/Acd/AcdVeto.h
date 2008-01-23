// $Header$
#ifndef CalibData_AcdVeto_h
#define CalibData_AcdVeto_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdVetoFitDesc
   *
   * @brief Description of a Veto threshold calibration.
   * 
   * This calibration consists of:
   *  - veto  = the Veto threshold in PHA counts (50% point of turn on curve)
   *  - width = the width of the turn on curve
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdVetoFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdVetoFitDesc& instance() {
      static const AcdVetoFitDesc desc;
      return desc;
    };     
  public:
    /// Trivial D'ctor
    virtual ~AcdVetoFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdVetoFitDesc()
      :AcdCalibDescription(AcdCalibData::VETO,"ACD_Veto"){
      addVarName("veto");
      addVarName("width");
    }
  };

  /** 
   * @class AcdVeto
   *
   * @brief A Veto threshold calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - veto  = the Veto threshold in PHA counts (50% point of turn on curve)
   *  - width = the width of the turn on curve
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdVeto : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_ThreshVeto;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::VETO;
    }
  public:
    /// Build from description and a set of values
    AcdVeto(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdVeto(float veto, float width, STATUS status) :
      AcdCalibObj(status,AcdVetoFitDesc::instance()){
      setVals(veto,width,status);
    }
    /// Trivial d'tor
    virtual ~AcdVeto() {}

    // Provide access to the values
    float getVeto() const { return (*this)[0];}
    float getWidth() const { return (*this)[1]; }
  };
}


#endif
