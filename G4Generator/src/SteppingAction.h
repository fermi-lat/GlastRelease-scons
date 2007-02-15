#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class IG4GeometrySvc;

class SteppingAction : public G4UserSteppingAction {

public:
    SteppingAction(IG4GeometrySvc* geoSvc);
    virtual ~SteppingAction(){};

    /// Action to be performed at each step
    virtual void UserSteppingAction(const G4Step*);

private:
    IG4GeometrySvc* m_geoSvc;
  
};
#endif
