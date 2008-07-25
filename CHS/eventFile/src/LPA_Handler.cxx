#include <stdio.h>

#include "eventFile/LPA_Handler.h"

namespace eventFile {

  // RSD method implementation
  unsigned int DgnHandlerRsdV0::gem_classes() const
  {
    return (status & DgnHandlerRsdV0::Gem_Classes_m) 
      >> DgnHandlerRsdV0::Gem_Classes_v;
  }

  unsigned int DgnHandlerRsdV0::stage() const
  {
    return status & DgnHandlerRsdV0::Stage_m;
  }

  unsigned int HipHandlerRsdV0::stage() const
  {
    return status & HipHandlerRsdV0::Stage_m;
  }

  unsigned int MipHandlerRsdV0::stage() const
  {
    return status & MipHandlerRsdV0::Stage_m;
  }

  // main implementation
  void LPA_Handler::dump( const char* pre, const char* post ) const
  {
    printf( "%stype      = %d (%s)%s", pre, type, typeName(), post );
    printf( "%sid        = %d (%s)%s", pre, id, handlerName(), post );
    printf( "%sversion   = %d%s", pre, version, post );
    printf( "%smasterKey = 0x%08X%s", pre, masterKey, post );
    printf( "%scfgKey    = 0x%08X%s", pre, cfgKey, post );
    printf( "%scfgId     = %d%s", pre, cfgId, post );
    printf( "%sstate     = 0x%08X (%s)%s", pre, state, stateName(), post );
    printf( "%sprescaler = %d (%s)%s", pre, prescaler, prescalerName(), post );
    const DgnHandlerRsdV0*      drsd( dgnRsdV0() );
    const GammaHandlerRsdV0*    grsd0( gammaRsdV0() );
    const GammaHandlerRsdV1*    grsd1( gammaRsdV1() );
    const GammaHandlerRsdV2*    grsd2( gammaRsdV2() );
    const HipHandlerRsdV0*      hrsd( hipRsdV0() );
    const MipHandlerRsdV0*      mrsd( mipRsdV0() );
    const PassthruHandlerRsdV0* prsd( passthruRsdV0() );
    if ( drsd ) {
      printf( "%s%s RSD:      status = 0x%08X%s", pre, handlerName(), drsd->status, post );
      printf( "%s%s RSD:      stage  = 0x%08X%s", pre, handlerName(), drsd->stage(), post );
    } else if ( grsd0 ) {
      printf( "%s%s RSD:    status  = 0x%08X%s", pre, handlerName(), grsd0->status, post );
      printf( "%s%s RSD:    stage   = 0x%08X%s", pre, handlerName(), grsd0->stage(), post );
      printf( "%s%s RSD:    energyValid  = %s%s", pre, handlerName(), ( (grsd0->energyValid) ? "True" : "False" ), post );
      if ( grsd0->energyValid ) {
	printf( "%s%s RSD:    energyInLeus = %d%s", pre, handlerName(), grsd0->energyInLeus, post );
      }
    } else if ( grsd1 ) {
      printf( "%s%s RSD:    status  = 0x%08X%s", pre, handlerName(), grsd1->status, post );
      printf( "%s%s RSD:    stage   = 0x%08X%s", pre, handlerName(), grsd1->stage(), post );
      printf( "%s%s RSD:    energyValid  = %s%s", pre, handlerName(), ( (grsd1->energyValid) ? "True" : "False" ), post );
      if ( grsd1->energyValid ) {
	printf( "%s%s RSD:    energyInLeus = %d%s", pre, handlerName(), grsd1->energyInLeus, post );
      }
    } else if ( grsd2 ) {
      printf( "%s%s RSD:    status  = 0x%08X%s", pre, handlerName(), grsd2->status, post );
      printf( "%s%s RSD:    stage   = 0x%08X%s", pre, handlerName(), grsd2->stage(), post );
      printf( "%s%s RSD:    energyValid  = %s%s", pre, handlerName(), ( (grsd2->energyValid) ? "True" : "False" ), post );
      if ( grsd2->energyValid ) {
	printf( "%s%s RSD:    energyInLeus = %d%s", pre, handlerName(), grsd2->energyInLeus, post );
      }
    } else if ( hrsd ) {
      printf( "%s%s RSD:      status = 0x%08X%s", pre, handlerName(), hrsd->status, post );
      printf( "%s%s RSD:      stage  = 0x%08X%s", pre, handlerName(), hrsd->stage(), post );
    } else if ( mrsd ) {
      printf( "%s%s RSD:      status = 0x%08X%s", pre, handlerName(), mrsd->status, post );
      printf( "%s%s RSD:      stage  = 0x%08X%s", pre, handlerName(), mrsd->stage(), post );
    } else if ( prsd ) {
      printf( "%s%s RSD:      status = 0x%08X%s", pre, handlerName(), prsd->status, post );
      printf( "%s%s RSD:      stage  = 0x%08X%s", pre, handlerName(), prsd->stage(), post );
    } else {
      printf( "%sno RSD for %s %s%s", pre, handlerName(), typeName(), post );
    }
  }

  const DgnHandlerRsdV0* LPA_Handler::dgnRsdV0() const
  {
    if ( has && version == 0 && id == LPA_Handler::DGN ) {
      return &(rsd.dgn0);
    }
    return NULL;
  }

  const GammaHandlerRsdV0* LPA_Handler::gammaRsdV0() const
  {
    if ( has && version == 0 && id == LPA_Handler::GAMMA ) {
      return &(rsd.gamma0);
    }
    return NULL;
  }

  const HipHandlerRsdV0* LPA_Handler::hipRsdV0() const
  {
    if ( has && version == 0 && id == LPA_Handler::HIP ) {
      return &(rsd.hip0);
    }
    return NULL;
  }

  const MipHandlerRsdV0* LPA_Handler::mipRsdV0() const
  {
    if ( has && version == 0 && id == LPA_Handler::MIP ) {
      return &(rsd.mip0);
    }
    return NULL;
  }

  const PassthruHandlerRsdV0* LPA_Handler::passthruRsdV0() const
  {
    if ( has && version == 0 && id == LPA_Handler::PASS_THRU ) {
      return &(rsd.passthru0);
    }
    return NULL;
  }

  const GammaHandlerRsdV1* LPA_Handler::gammaRsdV1() const
  {
    if ( has && version == 1 && id == LPA_Handler::GAMMA ) {
      return &(rsd.gamma1);
    }
    return NULL;
  }

  const GammaHandlerRsdV2* LPA_Handler::gammaRsdV2() const
  {
    if ( has && version == 2 && id == LPA_Handler::GAMMA ) {
      return &(rsd.gamma2);
    }
    return NULL;
  }

  static const char* HandlerTypeNames[] = { "Unknown",
					    "Filter",
					    "Monitor",
  };
  const char* LPA_Handler::typeName() const
  {
    int idx = type + abs( HandlerTypeMin );
    if ( idx >=0 && idx < HandlerTypeCnt ) {
      return HandlerTypeNames[ idx ];
    }
    return "Invalid";
  }

  static const char* HandlerIdNames[] = { "PassThru",
					  "Gamma",
					  "ACD Software Counters",
					  "Minimum Ionizing Particle",
					  "Highly Ionizing Particle",
					  "Diagnostic",
  };
  const char* LPA_Handler::handlerName() const
  {
    int idx = id + abs( HandlerIdMin );
    if ( idx >= 0 && idx < HandlerIdCnt ) {
      return HandlerIdNames[ id ];
    }
    return "Invalid";
  }

  static const char* RsdStateNames[] = { "Invalid",
					 "Passed",
					 "Suppressed",
					 "Vetoed",
					 "Leaked",
					 "Ignored",
  };
  const char* LPA_Handler::stateName() const
  {
    int idx = state + abs( RsdStateMin );
    if ( idx >= 0  && idx < RsdStateCnt ) {
      return RsdStateNames[ idx ];
    }
    return "Invalid";
  }

  static const char* LeakedPrescalerNames[] = { "Unsupported",
						"Input",
						"Output",
						"Cond00", "Cond01", "Cond02", "Cond03",
						"Cond04", "Cond05", "Cond06", "Cond07",
						"Cond08", "Cond09", "Cond10", "Cond11",
						"Cond12", "Cond13", "Cond14", "Cond15",
						"Cond16", "Cond17", "Cond18", "Cond19",
						"Cond20", "Cond21", "Cond22", "Cond23",
						"Cond24", "Cond25", "Cond26", "Cond27",
						"Cond28", "Cond29", "Cond30", "Cond31",
  };
  const char* LPA_Handler::prescalerName() const
  {
    int idx = prescaler + abs( LeakedPrescalerMin );
    if ( idx >= 0 && idx < LeakedPrescalerCnt ) {
      return LeakedPrescalerNames[ idx ];
    }
    return "Invalid";
  }
}
