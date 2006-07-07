#ifndef ENUMS_LSF_H
#define ENUMS_LSF_H

/*
 * The enums in this file are take from the DFI interface and should match with  
 * fsw/<version>/include/fsw/DFI/DfiEvent-*.hh
 *
 * Note that changes to this file affect the interpretation of persisted data.
 * So don't change it.
 *
*/

/// These are sentinel values that show that various LsfMetaEvent quantities are uninitialized
#define LSF_INVALID_ULONG 0xFFFFFFFFFFFFFFFF
#define LSF_INVALID_UINT 0xFFFFFFFF
#define LSF_INVALID_UCHAR 0xFF
#define LSF_INVALID_USHORT 0xFFFF
#define LSF_UNDEFINED 0xFFFF

namespace enums {

  namespace Lsf {
    
    namespace Open {
      /// The action that caused a datagram to be opened
      typedef enum { Start,        ///> First datagram after start of run
		     Resume,       ///> First datagram after a run was paused and then resumed
		     Continue,     ///> Any other datagram
		     Unspecified,  ///> Action information was lost
		     NumActions
      } Action;  

      /// The reason this datagram was opened
      typedef enum { Operator,     ///> Operator intervetion during a run
		     ModeChange,   ///> Automatic switching to or from TOO or ARR modes
		     TimedOut,     ///> Reached time limit requested by DAQ configuration
		     CountedOut,   ///> Reached event limit requested by DAQ configuration
		     Full,         ///> Datagram reached 64K size 
		     Unknown,      ///>
		     NumReasons    ///>
      } Reason;
    };
    
    namespace Close {
      
      /// The action that caused the datagram to be closed
      typedef enum { Stop,                 ///> Last datagram in a run which stopped normally
		     Abort,                ///> Last datagram in a run which was aborted
		     Pause,                ///> Last datagram before a run was paused
		     Continue,             ///> Any other datagram
		     Unspecified,          ///> Action information was lost
		     NumActions
      } Action;
      
      /// The reason this datagram was closed
      typedef enum { Operator,     ///> Operator intervetion during a run
		     ModeChange,   ///> Automatic switching to or from TOO or ARR modes
		     TimedOut,     ///> Reached time limit requested by DAQ configuration
		     CountedOut,   ///> Reached event limit requested by DAQ configuration
		     Full,         ///> Datagram reached 64K size 
		     Unknown,      ///>
		     NumReasons    ///>
      } Reason;
    };

    namespace TimeTone {
      enum { MISSING_GPS_MASK = 1,         ///> NO GPS lock, message w.r.t. LAT clock
	     MISSING_CPU_MASK = 2,         ///> NO 1-PPS signal at CPU level
	     MISSING_LAT_MASK = 4,         ///> NO 1-PPS signal at CPU level
	     MISSING_TIMETONE_MASK = 8,    ///> NO 1-PPS signal at Spacecraft 
	     EARLY_EVENT_MASK = 16         ///> Event arrived early
      };
    }
    
    /// The source of the data in this datagram
    typedef enum { NoCrate = -1, 
		   Epu0,
		   Epu1, 
		   Epu2, 
		   Siu0, 
		   Siu1, 
		   Aux,
		   Mixed, 
		   NumCrates
    } Crate;
    
    /// The operating mode that the data in the datagram was acquired under
    typedef enum { NoMode = -1, 
		   Normal,          ///> Standard operation
		   TOO,             ///> Target Of Oppurtunity
		   GRB0,            ///> GRB Suspected, but ARR not issued, not confirmed, or rejected
		   GRB1,            ///> ARR issued and accepted
		   GRB2,            ///> GRB is over, but Spacecraft still in ARR mode
		   Solar,           ///> Solar physics mode
		   Calibration,     ///> Charge injection for calibration
		   Diagnostic,      ///> Operator control for diagnostic
		   NumModes
    } Mode;
    
    /// The platform that a run was acquired on 
    typedef enum { NoPlatform = -1, 
		   Lat,         ///> The real thing
		   Testbed,     ///> hardware simulation
		   Host,        ///> software simulation
		   NumPlatforms 
    } Platform;
    
    /// The origin of this run
    typedef enum { NoOrigin = -1, 
		   Orbit,
		   MonteCarlo,
		   Ground,
		   NumOrigins
    } DataOrigin;
    
    /// The type of run.  
    typedef enum { NoRunType = -1, 
		   LPA,    ///> Particle data
		   AcdLCI, ///> Charge injection to Acd
		   CalLCI, ///> Charge injection to Cal
		   TkrLCI, ///> Charge injection to Tkr
		   NumRunTypes
    } RunType;

    typedef enum { AUTORANGE = 1, 
                   ZERO_SUPPRESSION = 2,
		   STROBE = 4,
    } Masks;

  }//namespace lsf
}//namespace enums


#endif
