// $Header$
#ifndef CalibData_AcdCheck_h
#define CalibData_AcdCheck_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdCheckDesc
   *
   * @brief Description of a check on various calibrations
   * 
   * This calibration consists of:
   *  - ped   = pedestal in mips should be 0
   *  - mip   = mip peak in mips should be 1
   *  - veto  = veto thresh in mips, should be 0.45
   *  - cno   = cno thresh in mips, should be 25
   *  - range = mismatch at range crossover, shoulbe be 0
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCheckDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdCheckDesc& instance() {
      static const AcdCheckDesc desc;
      return desc;
    };     
  public:
    /// Trivial D'ctor
    virtual ~AcdCheckDesc(){;};    
  private:    
    /// This is a singleton
    AcdCheckDesc()
      :AcdCalibDescription(AcdCalibData::MERITCALIB,"ACD_Check"){
      addVarName("ped");
      addVarName("mip");
      addVarName("veto");
      addVarName("cno");
      addVarName("range");
    }
  };

  /** 
   * @class AcdCheck
   *
   * @brief A calibration that summarizes how well other calibrations are working together
   * 
   * This calibration consists of:
   *  - ped   = pedestal in mips should be 0
   *  - mip   = mip peak in mips should be 1
   *  - veto  = veto thresh in mips, should be 0.45
   *  - cno   = cno thresh in mips, should be 25
   *  - range = mismatch at range crossover, shoulbe be 0
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCheck : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_Merit;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::MERITCALIB;
    }
  public:
    /// Build from description and a set of values
    AcdCheck(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdCheck(float ped, float mip, float veto, float cno, float range, STATUS status) :
      AcdCalibObj(status,AcdCheckDesc::instance()){
      setVals(ped,mip,veto,cno,range,status);
    }
    /// Trivial d'tor
    virtual ~AcdCheck() {}

    // Provide access to the values
    float getPed() const { return (*this)[0];}
    float getMip() const { return (*this)[1]; }
    float getVeto() const { return (*this)[2];}
    float getCno() const { return (*this)[3]; }
    float getRange() const { return (*this)[4]; }
  };
}


#endif
