#include "G4UIsession.hh"
#include "G4coutDestination.hh"
#include "globals.hh"

class UIsession : public G4UIsession
{
  // Base class of UI/GUI session
  
  public:
      virtual G4int ReceiveG4cout(G4String coutString){};
};
