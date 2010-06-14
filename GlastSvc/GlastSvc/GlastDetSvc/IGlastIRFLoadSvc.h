// $Id$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_IGlastIRFLoadSvc_
#define _H_IGlastIRFLoadSvc_


// includes
#include "Gaudi/Interfaces/IInterface.h"

// forward declarations
class   GlastDetector;

//!  This class defines the top-level interface to the GlastDetSvc class. 
//!  The concrete class will implement all of the methods underneath this
//! particular interface for the clients to use.
//
class   IGlastIRFLoadSvc : virtual public IInterface {
public:

    //! Try to load the next event from the IRF file. Returns
    //! a failure if the file is at EOF or if some other problem 
    //! arises.
    virtual StatusCode  loadNextEvent () = 0;

    //! Access the name of the IRF file that is being loaded.
    virtual std::string     irfFileName () = 0;

};

// IID declaration
extern const IID& IID_IGlastIRFLoadSvc;


#endif  // _H_IGlastDetSvc_
