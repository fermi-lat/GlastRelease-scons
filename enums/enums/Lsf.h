#ifndef ENUMS_LSF_H
#define ENUMS_LSF_H

namespace enums {
  namespace Lsf {

    /// The reason that a datagram was opened
    typedef enum { NoOpenReason = -1, 
		   Start=1, ///> first datagram afer start of run
		   Resume=2 ///> first datagram after a run was paused and then resumed
    } OpenReason;  

    /// The source of the request to open this datagram
    typedef enum { NoOpenRequest = -1, 
		   OpenOperator, 
		   OpenAutomatic, 
		   OpenUnknown 
    } OpenRequester;

    /// The source of the data in this datagram
    typedef enum { NoCrate = -1, 
		   Epu0, 
		   Epu1, 
		   Epu2, 
		   Siu0, 
		   Siu1, 
		   Aux 
    } Crate;

    /// The operating mode that the data in the datagram was acquired under
    typedef enum { NoMode = -1, 
		   Normal, 
		   TOO, ///> Target Of Oppurtunity
		   ARR, ///> Autonomous repointing
		   Calibration, 
		   Diagnostic 
    } Mode;
    
    /// The reason that a datagram was closed
    typedef enum { NoCloseReason = -1, 
		   Stop, 
		   Pause, 
		   Abort, 
		   Unknown 
    } CloseReason;

    /// The source of the request to open this datagram
    typedef enum { NoCloseRequest = -1, 
		   CloseOperator, 
		   CloseAutomatic, 
		   TimedOut, 
		   CountedOut, 
		   Full, 
		   CloseUnknown 
    } CloseRequester;
     
    /// The platform that a run was acquired on 
    typedef enum { NoPlatform = -1, 
		   Lat = 0, ///> The real thing
		   Testbed = 1, ///> hardware simulation
		   Host = 2 ///> software simulation
    } Platform;

    /// The origin of this run
    typedef enum { NoOrigin = -1, 
		   Orbit = 0, 
		   MonteCarlo = 1, 
		   Ground = 2 
    } DataOrigin;
 
    /// The type of run.  
    typedef enum { NoRunType = -1, 
		   LPA,    ///> Particle data
		   AcdLCI, ///> Charge injection to Acd
		   CalLCI, ///> Charge injection to Cal
		   TkrLCI  ///> Charge injection to Tkr
    } RunType;

  }//namespace lsf
}//namespace enums


#endif
