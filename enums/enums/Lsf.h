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

namespace enums {

  namespace Lsf {
    
    namespace Open {
      /// The action that caused a datagram to be opened
      typedef enum { Start,        ///> First datagram afer start of run
		     Resume,       ///> First datagram after a run was paused and then resumed
		     Continue,     ///> Any other datagram
		     Unspecified,  ///> Action information was lost
		     NumActions
      } Action;  

      /// The reason this datagram was opened
      typedef enum { Operator,     ///>
		     ModeChange,   ///>
		     TimedOut,     ///>
		     CountedOut,   ///>
		     Full,         ///>
		     Unknown,      ///>
		     NumReasons    ///>
      } Reason;
    };
    
    namespace Close {
      
      /// The action that caused the datagram to be closed
      typedef enum { Stop,                 ///> Last datagram in a run which stopped normally
		     Abort,                ///> Last datagram in a run which was aborted
		     Pause,                ///> Last datagram before a run was paused
		     Continue,             ///> 
		     Unspecified,          ///> Action information was lost
		     NumActions
      } Action;
      
      /// The reason this datagram was closed
      typedef enum { Operator,     ///>
		     ModeChange,   ///>
		     TimedOut,     ///>
		     CountedOut,   ///>
		     Full,         ///>
		     Unknown,      ///>
		     NumReasons    ///>
      } Reason;
    };
    
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
		   Normal, 
		   TOO,             ///> Target Of Oppurtunity
		   ARR,             ///> Autonomous repointing
		   Calibration, 
		   Diagnostic, 
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

  }//namespace lsf
}//namespace enums


#endif
