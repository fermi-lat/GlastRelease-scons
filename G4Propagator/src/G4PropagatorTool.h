// File and Version Information:
// $Header$
//

/** 
* @class G4PropagatorTool
*
* @brief A Gaudi Tool to provide an interface to the G4ParticlePropagator
*
* @author Tracy Usher
*
*/

#include "GlastSvc/Reco/IPropagatorTool.h"
#include "G4Generator/IG4GeometrySvc.h"
#include "G4ParticlePropagator.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

class G4PropagatorTool : public AlgTool, virtual public IPropagatorTool
{
 public: 

  G4PropagatorTool(const std::string& type, const std::string& name, const IInterface* parent);

  virtual ~G4PropagatorTool();
    
  /// Return a pointer to the G4ParticlePropagator singleton object
  virtual IKalmanParticle* getPropagator() {return G4ParticlePropagator::instance();}

  /// This is needed by propagator for getting at the geometry
  static G4VUserDetectorConstruction* UserDetector;

  /// This is needed for associating Geant4 volumes to Glast identifiers
  static IG4GeometrySvc::IdMap*       VolIdentMap;

  /// This is needed to access the Geant4 volume tracking
  static G4TransportationManager*     TransportationManager;
  
 private:
};
