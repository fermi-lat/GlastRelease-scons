#ifndef TackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class EventAction;

class TrackingAction : public G4UserTrackingAction {

  public:
    TrackingAction(EventAction* myEA);
    virtual ~TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
  
  private:
    EventAction* eventAction;
};

#endif
