// $Header$
#ifndef CalibData_CalibTest1_h
#define CalibData_CalibTest1_h

/** @class CalibTest1  

  Very simple calibration data-like class to be used for testing 
  calibration infrastructure

  @author J. Bogart
*/

#include "CalibData/CalibBase.h"
#include "CalibData/CalibModel.h"
// extern const CLID& CLID_Calib_CalibTest1;


namespace CalibData {
  class CalibTest1 : public CalibBase {
  public:
    CalibTest1(const std::string& name, int value, 
               const ITime& since, const ITime& till, int serNo = -1) :
      CalibBase(since, till, serNo), m_name(name), m_value(value) {}

    CalibTest1(const CalibTest1& other);

    virtual ~CalibTest1() {}

    // Re-implemented from DataObject
    /// Class ID of this instance
    inline virtual const CLID& clID() const { return classID(); } 
    
    /// Class ID of this class
    inline static  const CLID& classID() { return CLID_Calib_CalibTest1; };

    // Re-implemented from CalibBase
    virtual void    update(CalibTest1& other);
    std::string getValueName() const {return m_name;}
    int         getValue() const {return m_value;}

  private:
    std::string m_name;
    int         m_value;
  };

}


#endif
