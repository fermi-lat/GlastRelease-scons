
#ifndef CALDATA_H
#define CALDATA_H 1


#include "GaudiKernel/DataObject.h"

#include "data/CsIData.h"

/// Dummy class to make CsIData a DataObject as well

class CalData : virtual public DataObject , virtual public CsIData
{
    
};

        
#endif
