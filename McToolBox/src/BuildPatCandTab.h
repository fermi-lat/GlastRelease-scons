/**
 * @class BuildPatCandTab
 *
 * @brief This object builds the relational tables associating the Monte Carlo with the
 *        TkrRecon pattern recognition track candidates. 
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "GaudiKernel/IDataProviderSvc.h"

#ifndef BuildPatCandTab_h
#define BuildPatCandTab_h

namespace Event {

class BuildPatCandTab 
{
public:
    /// Standard Gaudi Tool interface constructor
    BuildPatCandTab(IDataProviderSvc* dataSvc);
   ~BuildPatCandTab();

private:
};

};

#endif
