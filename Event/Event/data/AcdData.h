
#ifndef ACDDATA_H
#define ACDDATA_H 1

#include "data/IVetoData.h"
#include "Gaudi/Kernel/DataObject.h"

/// Just to make IVetoData a DataObject

class AcdData : virtual public IVetoData , virtual public DataObject {

};

#endif
