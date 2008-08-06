// $Header$
#ifndef CalibData_AcdCnoFit_h
#define CalibData_AcdCnoFit_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdCnoFitFitDesc
   *
   * @brief Description of a CNO threshold calibration.
   * 
   * This calibration consists of:
   *  - cno   = the CNO threshold in PHA counts (50% point of turn on curve)
   *  - width = the width of the turn on curve
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCnoFitFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdCnoFitFitDesc& instance() {
      static const AcdCnoFitFitDesc desc;
      return desc;
    };        
  public:
    /// Trivial D'ctor
    virtual ~AcdCnoFitFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdCnoFitFitDesc()
      :AcdCalibDescription(AcdCalibData::CNO_FIT,"ACD_CnoFit"){
      addVarName("slope");
      addVarName("offset");
      addVarName("carbonPeak");
    }
  };

  /** 
   * @class AcdCnoFit
   *
   * @brief A CNO threshold calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - cno   = the CNO threshold in PHA counts (50% point of turn on curve)
   *  - width = the width of the turn on curve
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCnoFit : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_CnoFit;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::CNO_FIT;
    }
  public:
    /// Build from description and a set of values
    AcdCnoFit(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdCnoFit(float cno, float width, STATUS status) :
      AcdCalibObj(status,AcdCnoFitFitDesc::instance()){    
      setVals(cno,width,status);
    }
    /// Trivial d'tor
    virtual ~AcdCnoFit() {}

    // Provide access to the values
    float getSlope() const { return (*this)[0];}
    float getOffset() const { return (*this)[1]; }
    float getCarbonPeak() const { return (*this)[2]; }
  };
}


#endif
