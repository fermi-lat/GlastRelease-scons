// $Header$
#ifndef CalibData_AcdCno_h
#define CalibData_AcdCno_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdCnoFitDesc
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

  class AcdCnoFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdCnoFitDesc& instance() {
      static const AcdCnoFitDesc desc;
      return desc;
    };        
  public:
    /// Trivial D'ctor
    virtual ~AcdCnoFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdCnoFitDesc()
      :AcdCalibDescription(AcdCalibData::CNO,"ACD_Cno"){
      addVarName("cno");
      addVarName("width");
    }
  };

  /** 
   * @class AcdCno
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

  class AcdCno : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_ThreshHigh;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::CNO;
    }
  public:
    /// Build from description and a set of values
    AcdCno(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdCno(float cno, float width, STATUS status) :
      AcdCalibObj(status,AcdCnoFitDesc::instance()){    
      setVals(cno,width,status);
    }
    /// Trivial d'tor
    virtual ~AcdCno() {}

    // Provide access to the values
    float getCno() const { return (*this)[0];}
    float getWidth() const { return (*this)[1]; }
  };
}


#endif
