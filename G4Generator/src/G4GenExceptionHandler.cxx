//

#include "src/G4GenExceptionHandler.h"
#include "G4Generator/G4GenException.h"
#include "src/RunManager.h"
#include "G4StateManager.hh"
#include "G4ios.hh"
#include <stdlib.h>
#include "G4String.hh"

G4GenExceptionHandler::G4GenExceptionHandler() 
{
}

G4GenExceptionHandler::~G4GenExceptionHandler()
{
}

G4GenExceptionHandler::G4GenExceptionHandler(const G4GenExceptionHandler &)
:G4VExceptionHandler()
{
}

G4GenExceptionHandler& G4GenExceptionHandler::operator=(const G4GenExceptionHandler &)
{
   return *this;
}

G4int G4GenExceptionHandler::operator==(const G4GenExceptionHandler &right) const
{
   return (this == &right);
}

G4int G4GenExceptionHandler::operator!=(const G4GenExceptionHandler &right) const
{
   return (this != &right);
}

G4bool G4GenExceptionHandler::Notify(const char*         originOfException,
                                     const char*         exceptionCode,
                                     G4ExceptionSeverity severity,
                                     const char*         description)
{
    // Output initial error message
    G4cerr << G4endl;
    G4cerr << "*** G4Exception : " << exceptionCode << G4endl;
    G4cerr << "      issued by : " << originOfException << G4endl;
    G4cerr << description << G4endl;
  
    G4bool abortionForCoreDump = false;
    G4ApplicationState aps = G4StateManager::GetStateManager()->GetCurrentState();
  
    switch(severity)
    {
    case FatalException:
        G4cerr << "*** Fatal Exception *** core dump ***";
        abortionForCoreDump = true;
        break;
    case FatalErrorInArgument:
        G4cerr << "*** Fatal Error In Argument *** core dump ***";
        abortionForCoreDump = true;
        break;
    case RunMustBeAborted:
        if(aps==G4State_GeomClosed || aps==G4State_EventProc)
        {
            G4cerr << "*** Run Must Be Aborted ";
            RunManager::GetRunManager()->AbortRun(false);
        }
        abortionForCoreDump = false;
        break;
    case EventMustBeAborted:
        if(aps==G4State_EventProc)
        {
            std::string exceptionCodeStr(exceptionCode);

            G4cerr << "*** Event Must Be Aborted ";

            // Throw exception to our error handling so we can do book keeping
            // This requires a state reset first.
            G4StateManager::GetStateManager()->SetNewState(G4State_Idle);
            throw G4GenException("Aborting " + exceptionCodeStr + " Event");
        }

        abortionForCoreDump = false;
        break;
    default:
        G4cerr << "*** This is just a warning message.";
        abortionForCoreDump = false;
        break;
    }
  
    G4cerr << G4endl;
  
    return abortionForCoreDump;
}


