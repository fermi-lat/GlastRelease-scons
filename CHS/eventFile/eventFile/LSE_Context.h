/**
 * @class eventFile::LSE_Context
 *
 * @brief Class representing the acquisition context of the event.
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LSE_CONTEXT_HH
#define LSE_CONTEXT_HH

#include "eventFile/LSE_GemTime.h"

namespace eventFile {

  struct LSE_Context {

    struct FromCcsds {
      FromCcsds() {};
      void dump( const char* pre, const char* post ) const;
      int scid;
      int apid;
      double utc;
    };

    struct FromTimetone {
      FromTimetone() {};
      void dump( const char* pre, const char* post ) const;
      unsigned timeSecs;
      LSE_GemTime timeHack;
      unsigned incomplete;
      unsigned flywheeling;
      bool missingCpuPps;
      bool missingTimeTone;
      bool missingLatPps;
      bool missingGps;
    };

    struct FromScalers {
      FromScalers() {};
      void dump( const char* pre, const char* post ) const;
      unsigned long long elapsed;
      unsigned long long livetime;
      unsigned long long prescaled;
      unsigned long long discarded;
      unsigned long long sequence;
      unsigned long long deadzone;
    };

    struct FromRun {
      FromRun() {};
      void dump( const char* pre, const char* post ) const;
      int platform;
      int origin;
      unsigned groundId;
      unsigned startedAt;
      char platformTxt[16];
      char originTxt[16];
    };

    struct FromOpen {
      FromOpen() {};
      FromOpen( const DfiEvent::OpenInfo& open );
      void dump( const char* pre, const char* post ) const;
      unsigned modeChanges;
      unsigned datagrams;
      int action;
      int reason;
      int crate;
      int mode;
      char actionTxt[16];
      char reasonTxt[16];
      char crateTxt[16];
      char modeTxt[16];
    };

    struct FromClose {
      FromClose() {};
      FromClose( const DfiEvent::CloseInfo& close );
      void dump( const char* pre, const char* post ) const;
      int action;
      int reason;
      char actionTxt[16];
      char reasonTxt[16];
    };

    LSE_Context() {};
    LSE_Context( int scid, int apid, double utc, const DfiEvent::MetaEvent& meta );
    void dump() const;

    // from the enclosing CCSDS packet stream
    FromCcsds ccsds;

    // from the current timetone
    FromTimetone current;

    // from the previous timetone
    FromTimetone previous;

    // from the scalers
    FromScalers scalers;

    // from the run info
    FromRun run;

    // from the open info
    FromOpen open;

    // from the close info
    FromClose close;
  };

};

#endif
