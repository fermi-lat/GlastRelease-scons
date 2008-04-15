#include <stdio.h>

#include "eventFile/LPA_Handler.h"

namespace eventFile {

  void LPA_Handler::dump( const char* pre, const char* post ) const
  {
    printf( "%stype      = %d (%s)%s", pre, type, typeName(), post );
    printf( "%sid        = %d (%s)%s", pre, id, handlerName(), post );
    printf( "%smasterKey = 0x%08X%s", pre, masterKey, post );
    printf( "%scfgKey    = 0x%08X%s", pre, cfgKey, post );
    printf( "%scfgId     = %d%s", pre, cfgId, post );
    printf( "%sstate     = 0x%08X%s", pre, state, post );
    printf( "%sprescaler = %d%s", pre, prescaler, post );
    const DgnHandlerRsdV0*      drsd( dgnRsdV0() );
    const GammaHandlerRsdV0*    grsd( gammaRsdV0() );
    const HipHandlerRsdV0*      hrsd( hipRsdV0() );
    const MipHandlerRsdV0*      mrsd( mipRsdV0() );
    const PassthruHandlerRsdV0* prsd( passthruRsdV0() );
    if ( drsd ) {
      printf( "%s%s RSD:      status = 0x%08X%s", pre, handlerName(), drsd->status, post );
    } else if ( grsd ) {
      printf( "%s%s RSD:    status  = 0x%08X%s", pre, handlerName(), grsd->status, post );
      printf( "%s%s RSD:    stage   = 0x%08X%s", pre, handlerName(), grsd->stage, post );
      printf( "%s%s RSD:    energyValid  = %s%s", pre, handlerName(), ( (grsd->energyValid) ? "True" : "False" ), post );
      if ( grsd->energyValid ) {
	printf( "%s%s RSD:    energyInLeus = %d%s", pre, handlerName(), grsd->energyInLeus, post );
      }
    } else if ( hrsd ) {
      printf( "%s%s RSD:      status = 0x%08X%s", pre, handlerName(), hrsd->status, post );
    } else if ( mrsd ) {
      printf( "%s%s RSD:      status = 0x%08X%s", pre, handlerName(), mrsd->status, post );
    } else if ( prsd ) {
      printf( "%s%s RSD: status = 0x%08X%s", pre, handlerName(), prsd->status, post );
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

  static const char* HandlerTypeNames[] = { "Unknown",
					    "Filter",
					    "Monitor",
  };
  const char* LPA_Handler::typeName() const
  {
    if ( type+1 < HandlerTypeCnt ) {
      return HandlerTypeNames[ type+1 ];
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
    if ( id >= 0 && id < MaxHandlerIds ) {
      return HandlerIdNames[ id ];
    }
    return "Invalid";
  }

}
