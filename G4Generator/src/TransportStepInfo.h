#ifndef TransportStepInfo_h
#define TransportStepInfo_h 1

#include "G4VPhysicalVolume.hh"
#include "globals.hh"
#include <vector>

/** 
* @class TransportStepInfo
*
* @brief Class internal to ParticleTransporter for keeping track of an
* individual step
*
* @author Tracy Usher
*
*/

//Internal class for keeping track of step info
class TransportStepInfo
{
public:
  TransportStepInfo(const G4ThreeVector& coords, G4double step, 
                    G4VPhysicalVolume* pVolume) :
      position(coords),arcLen(step),pCurVolume(pVolume) {};
   ~TransportStepInfo() {};

    G4ThreeVector       GetCoords() {return position;}
    G4double            GetArcLen() {return arcLen;}
    G4VPhysicalVolume*  GetVolume() {return pCurVolume;}

    void SetCoords(const G4ThreeVector& coords) {position   = coords;}
    void SetArcLen(const G4double newArcLen)    {arcLen     = newArcLen;}
    void                SetVolume(G4VPhysicalVolume* volume)   {pCurVolume = volume;}
    
private:
    G4ThreeVector      position;
    G4double           arcLen;
    G4VPhysicalVolume* pCurVolume;
};

typedef std::vector<TransportStepInfo*>                 StepVector;
typedef std::vector<TransportStepInfo*>::iterator       StepPtr;
typedef std::vector<TransportStepInfo*>::const_iterator ConstStepPtr;

#endif




