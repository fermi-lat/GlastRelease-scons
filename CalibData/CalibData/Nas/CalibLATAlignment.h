// $Header$
#ifndef CalibData_LATAlignment_h
#define CalibData_LATAlignment_h

/** @class CalibTest1  

  Very simple calibration data-like class to be used for testing 
  calibration infrastructure

  @author J. Bogart
*/

#include "CalibData/CalibBase.h"
#include "CalibData/CalibModel.h"
#include <vector>

// extern const CLID& CLID_Calib_NAS_LATAlignment;


namespace CalibData {
  typedef double ALIGN_ROT[3];
  class CalibLATAlignment : public CalibBase {
  public:
    CalibLATAlignment(double rx, double ry, double rz, 
                      const std::string& units,
                      const ITime& since, const ITime& till,int serNo);

    CalibLATAlignment(const CalibLATAlignment& other);

    virtual ~CalibLATAlignment() {}

    // Re-implemented from DataObject
    inline virtual const CLID& clID() const { return classID(); } 
    
    inline static  const CLID& classID() { return CLID_Calib_NAS_LATAlignment; };


    // Re-implemented from CalibBase
    virtual StatusCode   update(CalibBase& other, MsgStream* log);

    double getRx() const {return m_r[0];}; //Rotation around x-axis
    double getRy() const {return m_r[1];}; //Rotation around y-axis
    double getRz() const {return m_r[2];}; //Rotation around z-axis

    const ALIGN_ROT* getR() {
      return &m_r;
    }
    std::string getUnits() {return m_units;}
    
    
  protected:

  private:
    //    double m_r[3];
    ALIGN_ROT m_r;
    std::string m_units;
  };
}


#endif
