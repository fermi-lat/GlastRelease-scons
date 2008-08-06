// $Header$
#ifndef CalibData_AcdCarbon_h
#define CalibData_AcdCarbon_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdCarbonFitDesc
   *
   * @brief Description of Carbon (aka MIP peak) calibration
   * 
   * This calibration consists of:
   *  - peak  = the mip peak in PHA counts above pedestal
   *  - width = the width of the MIP peak in PHA counts
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCarbonFitDesc : public AcdCalibDescription {
  public:
    /// Get this description
   static const AcdCarbonFitDesc & instance() {
      static const AcdCarbonFitDesc desc;
      return desc;
    }
  public:
    /// Trivial D'ctor
    virtual ~AcdCarbonFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdCarbonFitDesc()
      :AcdCalibDescription(AcdCalibData::CARBON,"ACD_Carbon"){
      addVarName("peak");
      addVarName("width");
    }
  };

  /** 
   * @class AcdCarbon
   *
   * @brief A Carbon (aka MIP peak) calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - peak  = the mip peak in PHA counts above pedestal
   *  - width = the width of the MIP peak in PHA counts
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCarbon : public AcdCalibObj {    
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_Carbon;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::GAIN;
    }
  public:
    /// Build from description and a set of values
    AcdCarbon(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdCarbon(float peak, float width, STATUS status) :
      AcdCalibObj(status,AcdCarbonFitDesc::instance()){
      setVals(peak,width,status);
    }
     /// Trivial d'tor
    virtual ~AcdCarbon() {}

    // Provide access to the values
    float getPeak() const { return (*this)[0];}
    float getWidth() const { return (*this)[1]; }
  };
}

#endif
