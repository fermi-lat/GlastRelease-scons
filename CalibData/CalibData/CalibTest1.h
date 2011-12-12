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


namespace CalibData {
  class CalibTest1 : public CalibBase {
  public:
    CalibTest1(const std::string& name, int value, 
               const Gaudi::Time& since, const Gaudi::Time& till, int serNo = -1);

    virtual ~CalibTest1() {}

    // Re-implemented from DataObject
    inline virtual const CLID& clID() const { return classID(); } 
    
    inline static  const CLID& classID() { return CLID_Calib_CalibTest1; };


    // Re-implemented from CalibBase
    virtual StatusCode   update(CalibBase& other, MsgStream* log);

    std::string getValueName() const;
    int         getValue() const {return m_value;}

  protected:

  private:
    std::string m_name;
    int         m_value;
  };
}


#endif
