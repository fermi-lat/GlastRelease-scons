/**
 * @class McBuildPatCandTab
 *
 * @brief Represents a Monte Carlo hit in a tracker layer
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ContainedObject.h"

#ifndef McBuildPatCandTab_h
#define McBuildPatCandTab_h

namespace Event {

class McBuildPatCandTab : virtual public ContainedObject
{
public:
    /// Standard Gaudi Tool interface constructor
    McBuildPatCandTab(DataSvc* dataSvc);
   ~McBuildPatCandTab();

private:
};

};

#endif