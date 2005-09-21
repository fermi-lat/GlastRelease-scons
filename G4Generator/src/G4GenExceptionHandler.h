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

#ifndef G4GenExceptionHandler_h
#define G4GenExceptionHandler_h 1

#include "globals.hh"
#include "G4VExceptionHandler.hh"
#include "G4ExceptionSeverity.hh"

class G4GenExceptionHandler : public G4VExceptionHandler
{

public:

  G4GenExceptionHandler();
  virtual ~G4GenExceptionHandler();
  G4int operator==(const G4GenExceptionHandler &right) const;
  G4int operator!=(const G4GenExceptionHandler &right) const;

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

  G4GenExceptionHandler(const G4GenExceptionHandler &right);
  G4GenExceptionHandler& operator=(const G4GenExceptionHandler &right);

};

#endif
