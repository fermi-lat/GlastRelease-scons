// File and Version Information:
// $Header$
//
// Description: Implementation of a G4Exception handler. Intermediate between G4 and 
//              G4Generator exception handling. 
//      
// Author(s):
//      T. Usher
//
// ------------------------------------------------------------

#ifndef G4PropagatorExceptionHandler_h
#define G4PropagatorExceptionHandler_h 1

#include "globals.hh"
#include "G4VExceptionHandler.hh"
#include "G4ExceptionSeverity.hh"

class G4PropagatorExceptionHandler : public G4VExceptionHandler
{

public:

  G4PropagatorExceptionHandler();
  virtual ~G4PropagatorExceptionHandler();
  G4int operator==(const G4PropagatorExceptionHandler &right) const;
  G4int operator!=(const G4PropagatorExceptionHandler &right) const;

public: // with description

  virtual G4bool Notify(const char* originOfException,
                        const char* exceptionCode,
                        G4ExceptionSeverity severity,
                        const char* description);
    // Virtual method which will be invoked by G4StateManager when
    // G4Exception occurs.
    // If TRUE returned, core dump will be generated, while FALSE returned,
    // program execution continues.

private:

  G4PropagatorExceptionHandler(const G4PropagatorExceptionHandler &right);
  G4PropagatorExceptionHandler& operator=(const G4PropagatorExceptionHandler &right);

};

#endif
