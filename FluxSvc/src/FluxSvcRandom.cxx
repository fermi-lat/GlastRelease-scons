
#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"
#include "GlastSvc/GlastRandomSvc/RandomAccess.h"
#include "GaudiKernel/ToolFactory.h"

class FluxSvcRandom : public RandomAccess {
public:
  FluxSvcRandom( const std::string& type, const std::string& name, const
		  IInterface* parent):RandomAccess(type,name,parent) {}
};

// Static factory for instantiation of algtool objects
static ToolFactory<FluxSvcRandom> s_factory;
const IToolFactory& FluxSvcRandomFactory = s_factory;
