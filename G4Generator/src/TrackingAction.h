#ifndef TackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

/** 
 * @class TrackingAction
 *
 * @brief A class to help the filling of the McParticle tree
 *
 * This class uses standard Geant4 mechanisms, derived from
 * G4UserTrackingAction, to intercept creation of new tracks to fill the
 * McParticle tree
 *  
 * @author R.Giannitrapani
 *
 * $Header$
 */
class TrackingAction : public G4UserTrackingAction {

  public:
    TrackingAction();
    virtual ~TrackingAction(){};

    /// Action to be performed before a new track is created
    virtual void PreUserTrackingAction(const G4Track*);
    /// Action to be performed after a track is dead
    virtual void PostUserTrackingAction(const G4Track*);
  
  private:
};
#endif
