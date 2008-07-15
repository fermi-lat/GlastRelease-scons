// $Header$
#ifndef CalibData_AcdRibbon_h
#define CalibData_AcdRibbon_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdRibbonFitDesc
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

  class AcdRibbonFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdRibbonFitDesc& instance() {
      static const AcdRibbonFitDesc desc;
      return desc;
    };        
  public:
    /// Trivial D'ctor
    virtual ~AcdRibbonFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdRibbonFitDesc()
      :AcdCalibDescription(AcdCalibData::RIBBON,"ACD_Ribbon"){
      addVarName("Frac_n3");
      addVarName("Frac_n2");
      addVarName("Frac_n1");
      addVarName("Frac_p1");
      addVarName("Frac_p2");
      addVarName("Frac_p3");      
      addVarName("Norm");
    }
  };

  /** 
   * @class AcdRibbon
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

  class AcdRibbon : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_Ribbon;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::RIBBON;
    }
  public:
    /// Build from description and a set of values
    AcdRibbon(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdRibbon(float n3, float n2, float n1, float p1, float p2, float p3, float norm, STATUS status) :
      AcdCalibObj(status,AcdRibbonFitDesc::instance()){    
      setVals(n3,n2,n1,p1,p2,p3,norm,status);
    }
    /// Trivial d'tor
    virtual ~AcdRibbon() {}

    // Provide access to the values
    float getN1() const { return (*this)[0];}
    float getN2() const { return (*this)[1];}
    float getN3() const { return (*this)[2];}
    float getP1() const { return (*this)[3];}
    float getP2() const { return (*this)[4];}
    float getP3() const { return (*this)[5];}
    float getNorm() const { return (*this)[6];}
  };
}


#endif
