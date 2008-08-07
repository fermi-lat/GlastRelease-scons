// $Header$
#ifndef CalibData_AcdPE_h
#define CalibData_AcdPE_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdPEFitDesc
   *
   * @brief Description of the number of photoelectrons/mip calibration
   *  - PE   = the number of photoelectrons/ mip
   * 
   * This calibration consists of:
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdPEFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdPEFitDesc& instance() {
      static const AcdPEFitDesc desc;
      return desc;
    };        
  public:
    /// Trivial D'ctor
    virtual ~AcdPEFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdPEFitDesc()
      :AcdCalibDescription(AcdCalibData::PE,"ACD_PE"){
      addVarName("PE");
    }
  };

  /** 
   * @class AcdPE
   *
   * @brief Number of photoelectrons/mip calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - PE   = the number of photoelectrons/ mip
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdPE : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_PE;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::PE;
    }
  public:
    /// Build from description and a set of values
    AcdPE(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdPE(float PE, STATUS status) :
      AcdCalibObj(status,AcdPEFitDesc::instance()){    
      setVals(PE,status);
    }
    /// Trivial d'tor
    virtual ~AcdPE() {}

    // Provide access to the values
    float getPE() const { return (*this)[0];}
  };
}


#endif
