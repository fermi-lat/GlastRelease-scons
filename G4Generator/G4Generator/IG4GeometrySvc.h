#ifndef IG4GEOMETRYSVC_H
#define IG4GEOMETRYSVC_H
// File and Version Information:
// $Header$
//
// Description: Defines interface for Geant4 Geometry Service
//
// Author(s):
//      T.Usher

#include "GaudiKernel/IInterface.h"
#include <map>

class G4VUserDetectorConstruction;
class G4TransportationManager;
class G4VPhysicalVolume;
#include "idents/VolumeIdentifier.h"

/** 
 * @class IG4GeometrySvc
 *
 * @brief A Gaudi Service for setting up the geometry for Geant4 and the Propagator
 *
 * @author Tracy Usher
 *
 */

// Declaration of the interface ID
static const InterfaceID IID_IG4GeometrySvc("IG4GeometrySvc", 1 , 0); 

class IG4GeometrySvc : virtual public IInterface
{
 public: 
  /// Return pointer to the constructed Detector
  virtual G4VUserDetectorConstruction* getDetector()= 0;

  /// Return pointer to the G4 Transportation Manager
  virtual G4TransportationManager* getTransportationManager() = 0;

  /// Return pointer to the ID map
  typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier > IdMap;
  virtual IdMap* getIdMap() = 0;

  /// queryInterface - for implementing a Service this is necessary
  static const InterfaceID& interfaceID() { return IID_IG4GeometrySvc; }
};

#endif