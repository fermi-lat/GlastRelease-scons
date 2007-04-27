#include <stdio.h>

#include "eventFile/LSE_GemTime.h"
#include "eventFile/LSE_Context.h"

namespace eventFile {

  void FromCcsds::dump( const char* pre, const char* post ) const 
  {
    printf( "%sscid = %d%s", pre, scid, post );
    printf( "%sapid = %d%s", pre, apid, post );
    printf( "%sutc  = %18.6f%s", pre, utc, post );
  }

  void FromTimetone::dump( const char* pre, const char* post ) const
  {
    printf( "%ssecs            = 0x%08x (%d)%s", pre, timeSecs, timeSecs, post );
    printf( "%sincomplete      = 0x%08x%s", pre, incomplete, post );
    printf( "%sflywheeling     = 0x%08x%s", pre, flywheeling, post );
    printf( "%smissingTimeTone = %s%s", pre, (missingTimeTone) ? "True" : "False", post );
    printf( "%smissingCpuPps   = %s%s", pre, (missingCpuPps)   ? "True" : "False", post );
    printf( "%smissingLatPps   = %s%s", pre, (missingLatPps)   ? "True" : "False", post );
    printf( "%searlyEvent      = %s%s", pre, (earlyEvent)      ? "True" : "False", post );
    printf( "%smissingGps      = %s%s", pre, (missingGps)      ? "True" : "False", post );
    timeHack.dump( "   ", "\n" );
  }

  void FromScalers::dump( const char* pre, const char* post ) const
  {
    printf( "%selapsed   = 0x%016llx = %llu%s", pre, elapsed, elapsed, post );
    printf( "%slivetime  = 0x%016llx = %llu%s", pre, livetime, livetime, post );
    printf( "%sprescaled = 0x%016llx = %llu%s", pre, prescaled, prescaled, post );
    printf( "%sdiscarded = 0x%016llx = %llu%s", pre, discarded, discarded, post );
    printf( "%ssequence  = 0x%016llx = %llu%s", pre, sequence, sequence, post );
    printf( "%sdeadzone  = 0x%016llx = %llu%s", pre, deadzone, deadzone, post );
  }

  void FromRun::dump( const char* pre, const char* post ) const
  {
    printf( "%sgroundid = 0x%06X, started = 0x%08X (%d )%s", pre, groundId, startedAt, startedAt, post );
    printf( "%splatform = %s (%d)%s", pre, platformTxt, platform, post );
    printf( "%sorigin   = %s (%d)%s", pre, originTxt, origin, post );
  }

  void FromOpen::dump( const char* pre, const char* post ) const
  {
    printf( "%snmodes = %d, ndgms = %d%s", pre, modeChanges, datagrams, post );
    printf( "%saction = %s (%d)%s", pre, actionTxt, action, post );
    printf( "%sreason = %s (%d)%s", pre, reasonTxt, reason, post );
    printf( "%scrate  = %s (%d)%s", pre, crateTxt, crate, post );
    printf( "%smode   = %s (%d)%s", pre, modeTxt, mode, post );
  }

  void FromClose::dump( const char* pre, const char* post ) const
  {
    printf( "%saction = %s (%d)%s", pre, actionTxt, action, post );
    printf( "%sreason = %s (%d)%s", pre, reasonTxt, reason, post );
  }

  void LSE_Context::dump() const
  {
    ccsds.dump(    " ccsds:    ", "\n" );
    current.dump(  " current:  ", "\n" );
    previous.dump( " previous: ", "\n" );
    scalers.dump(  " scalers:  ", "\n" );
    run.dump(      " run:      ", "\n" );
    open.dump(     " open:     ", "\n" );
    close.dump(    " close:    ", "\n" );
  }

}
