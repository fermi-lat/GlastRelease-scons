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
#ifdef WIN32
#define LSF_INVALID_ULONG 0xFFFFFFFFFFFFFFFF
#else
#define LSF_INVALID_ULONG 0xFFFFFFFFFFFFFFFFULL
#endif
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
      enum { SOURCE_GPS_MASK = 1,         ///> GPS lock, message w.r.t. LAT clock
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

    // The type of translated keys.
    typedef enum { NoKeysType = -1,
		   LpaKeys,   ///> master, ignore, list-of-LPA_DB
		   LciKeys,   ///> master, ignore, script
		   NumKeysTypes
    } KeysType;

    // CAL readout bits
    typedef enum { AUTORANGE = 1, 
                   ZERO_SUPPRESSION = 2,
		   STROBE = 4,
    } Masks;


       /** enumerate the types of handlers */
    typedef enum {
      Unknown = -1,
      Filter  =  0,
      Monitor =  1,

      // management values
      HandlerTypeMin = -1,
      HandlerTypeMax =  1,
      HandlerTypeCnt = HandlerTypeMax - HandlerTypeMin + 1,
    } HandlerType;

    /** enumerate the handler identities */
    typedef enum {
      PASS_THRU = 0,
      GAMMA     = 1,
      ASC       = 2,
      MIP       = 3,
      HIP       = 4,
      DGN       = 5,
      MaxHandlerIds, // deprecated
      
      // Management values
      HandlerIdMin=0,
      HandlerIdMax=5,
      HandlerIdCnt=HandlerIdMax-HandlerIdMin+1,
    } HandlerId;

    /** enumerate the overall result states for an event */
    typedef enum {
      INVALID    = -1,  /// Invalid state
      PASSED     =  0,  /// Event passed the filter
      SUPPRESSED =  1,  /// Event passed but was suppressed by an output prescaler
      VETOED     =  2,  /// Event vetoed by the filter
      LEAKED     =  3,  /// Event vetoed but leaked by an output prescaler
      IGNORED    =  4,  /// Event ignored due to input prescaler
   
      // Management values
      RsdStateMin = -1,
      RsdStateMax = 4,
      RsdStateCnt = RsdStateMax-RsdStateMin+1,
    } RsdState;

    /** enumerate valid values for the prescaler data member 
	note the values > 0 indicate which condition-specifc
	prescaler expired.
    */
    typedef enum {
      UNSUPPORTED = -3, /// prescaler info not available (always for RSD V0)
      INPUT       = -2, /// event not analyzed due to input prescaler
      OUTPUT      = -1, /// event leaked by overall output prescaler

      /** event leaked by condition-specific prescaler */
      COND00, COND01, COND02, COND03, COND04, COND05, COND06, COND07,
      COND08, COND09, COND10, COND11, COND12, COND13, COND14, COND15,
      COND16, COND17, COND18, COND19, COND20, COND21, COND22, COND23,
      COND24, COND25, COND26, COND27, COND28, COND29, COND30, COND31,

      // management values
     LeakedPrescalerMin = -3,
     LeakedPrescalerMax = 31,
     LeakedPrescalerCnt = LeakedPrescalerMax - LeakedPrescalerMin + 1, 
    } LeakedPrescaler;




  }//namespace lsf
}//namespace enums


#endif
