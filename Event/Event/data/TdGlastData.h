// TdGlastData.h: interface for the TdGlastData class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TdGlastData_H
#define TdGlastData_H 1

#include "data/SiData.h"
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include "geometry/Point.h"
#include "idents/ModuleId.h"

#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "data/GlastData.h"


// forward declarations

#include <iostream>
#include <vector>

extern const CLID& CLID_TdGlastData;

/*! \class TdGlastData
    
      This class is present simply to introduce DataObject into the LdGlastData.
      This modifocation was done mainly for the purpose of aoRecon however it
      also simplifies the loading of TdSiData and TdCsIData objects.

*/


class TdGlastData : virtual public GlastData , virtual public DataObject {

public:


};

#endif
