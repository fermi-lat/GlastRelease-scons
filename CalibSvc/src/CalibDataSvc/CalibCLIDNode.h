// $Header$

#ifndef CalibCLIDNode_H
#define CalibCLIDNode_H

#include <iostream>
#include <string>
//#include "CalibCnv/ICalibCnvSvc.h"

/** @class CalibCLIDNode

  Trivial DataObject.  Only extra data is a field to contain class id
  of child nodes (which will contain actual calibration data set).

@author Joanne Bogart
 $Header$

*/

// extern const CLID& CLID_Calib_CalibCLIDNode;

class CalibCLIDNode : virtual public DataObject {

  public:

    CalibCLIDNode(const CLID childClassID) :
      DataObject(), m_childClassID(childClassID) {}

    // Having these inline in include file could cause problems, in
    // fact the static member could already be a problem.  Will code
    // linked into different shareables have different copies of
    // 
    virtual const CLID& clID() const {return CalibCLIDNode::classID(); }
    static const CLID& classID();

    inline CLID getChildClassID() const {return m_childClassID;}

    virtual std::ostream& fillStream(std::ostream& s) const;

  private:
    CLID  m_childClassID;
  };
}





