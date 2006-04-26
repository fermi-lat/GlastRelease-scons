#ifndef TholdCIMgr_H
#define TholdCIMgr_H
// $Header$
// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalTholdCI.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB
// STD

using namespace CalUtil;
using namespace idents;

using namespace CalibData;

class CalCalibSvc;

/** @class TholdCIMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal charge-injection measured threshold data.
*/

class TholdCIMgr : public CalibItemMgr {
 public:
  TholdCIMgr(CalCalibShared &ccsShared) 
    : CalibItemMgr(CAL_TholdCI, ccsShared)
    {};

  const CalTholdCI *getTholdCI(FaceIdx faceIdx);

 private:
  StatusCode loadIdealVals();

  
  StatusCode genLocalStore();
  
  /// Validate TDS data entry (for empty ptrs & fun stuff like that)
  bool validateRangeBase(CalTholdCI *tholdCI);

  auto_ptr<CalTholdCI> m_idealTholdCI;
  
};

#endif
