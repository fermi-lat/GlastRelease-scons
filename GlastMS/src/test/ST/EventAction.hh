// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

// --------------------------------------------------------------
//      GEANT 4 - GLAST  2000
//---------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"  
#include "globals.hh"

//class EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetPrintModulo(G4int    val)  {printModulo = val;};
  
    void SetStartEn(G4double val) {startEn = val;};

    void SetEndEn(G4double val) {endEn = val;};

    void SetInitialDir(const G4ThreeVector& aValue)
      {initialDir = aValue;};
    void SetFinalDir(const G4ThreeVector& aValue)
      {finalDir = aValue;};

  private:
    G4int       layerCollID;                
    G4int       printModulo;                         
    G4double	startEn;
    G4double    endEn;
    
    G4ThreeVector initialDir;
    G4ThreeVector finalDir;
};

#endif

    
