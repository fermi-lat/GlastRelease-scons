// Implementation file for RandomAccess which sets and gets 
// random engines
//
// $Header
//
// Author: Toby Burnett, Karl Young


#include "GlastSvc/GlastRandomSvc/RandomAccess.h"
#include "GaudiKernel/ToolFactory.h"


static ToolFactory<RandomAccess> s_factory;
const IToolFactory& RandomAccessFactory = s_factory;

