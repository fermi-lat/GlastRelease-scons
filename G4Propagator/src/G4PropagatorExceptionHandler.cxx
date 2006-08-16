//

#include "G4PropagatorExceptionHandler.h"
#include "G4Generator/G4GenException.h"
#include "src/RunManager.h"
#include "G4StateManager.hh"
#include "G4ios.hh"
#include <stdlib.h>
#include "G4String.hh"

G4PropagatorExceptionHandler::G4PropagatorExceptionHandler() 
{
}

G4PropagatorExceptionHandler::~G4PropagatorExceptionHandler()
{
}

G4PropagatorExceptionHandler::G4PropagatorExceptionHandler(const G4PropagatorExceptionHandler &)
:G4VExceptionHandler()
{
}

G4PropagatorExceptionHandler& G4PropagatorExceptionHandler::operator=(const G4PropagatorExceptionHandler &)
{
   return *this;
}

G4int G4PropagatorExceptionHandler::operator==(const G4PropagatorExceptionHandler &right) const
{
   return (this == &right);
}

G4int G4PropagatorExceptionHandler::operator!=(const G4PropagatorExceptionHandler &right) const
{
   return (this != &right);
}

G4bool G4PropagatorExceptionHandler::Notify(const char*  originOfException,
                                     const char*         exceptionCode,
                                     G4ExceptionSeverity severity,
                                     const char*         description)
{
    // Output initial error message
    G4cerr << G4endl;
    G4cerr << "*** G4Exception : " << exceptionCode << G4endl;
    G4cerr << "      issued by : " << originOfException << G4endl;
    G4cerr << description << G4endl;

    std::string exceptString = "G4Exception: " + std::string(exceptionCode) + "\n"
                             + "issued by:   " + std::string(originOfException) ;

    throw G4GenException("Handling bad event: \n" + exceptString);
  
    G4bool abortionForCoreDump = false;
  
    G4cerr << G4endl;
  
    return abortionForCoreDump;
}


