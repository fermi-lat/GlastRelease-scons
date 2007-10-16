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
  class CalibLATAlignment : public CalibBase {
  public:
    CalibLATAlignment(double phi, double theta, double psi,
                      const ITime& since, const ITime& till,int serNo);

    CalibLATAlignment(const CalibLATAlignment& other);

    virtual ~CalibLATAlignment() {}

    // Re-implemented from DataObject
    inline virtual const CLID& clID() const { return classID(); } 
    
    inline static  const CLID& classID() { return CLID_Calib_NAS_LATAlignment; };


    // Re-implemented from CalibBase
    virtual StatusCode   update(CalibBase& other, MsgStream* log);

    double getRoll() const {return m_roll;}; //Rotation around x-axis
    double getPitch() const {return m_pitch;}; //Rotation around y-axis
    double getYaw() const {return m_yaw;}; //Rotation around z-axis
    
  protected:

  private:
    double m_roll;
    double m_pitch;
    double m_yaw;
  };
}


#endif
