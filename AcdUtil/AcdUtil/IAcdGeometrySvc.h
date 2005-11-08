/** @file IAcdGeometrySvc.h
 @brief Abstract interface to TkrGeometrySvc (q.v.)

  $Header$
*/

#ifndef __IACDGEOMETRYSVC_H
#define __IACDGEOMETRYSVC_H 1

#include "GaudiKernel/IInterface.h"


#include "CLHEP/Geometry/Point3D.h"

#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "geometry/Ray.h"

#include "AcdUtil/AcdDetectorList.h"

#include <string>

/** 
 * @class IAcdGeometrySvc
 *
 * @brief Abstract interface to AcdGeometrySvc 
 * 
 * @author Heather Kelly 
 */

static const InterfaceID IID_IAcdGeometrySvc("IAcdGeometrySvc", 1, 1); 

class IAcdGeometrySvc : public virtual IInterface
{
public:

    //! Constructor of this form must be provided

    static const InterfaceID& interfaceID() { return IID_IAcdGeometrySvc; }
    
    //Retrieve stored information

    virtual int    numTiles()     const = 0;
    virtual int    numRibbons()   const = 0;


    virtual StatusCode getDimensions(const idents::VolumeIdentifier &volIId,
                 std::vector<double> &dims, HepPoint3D &xT) const = 0;

    virtual StatusCode getDetectorDimensions(
                 const idents::VolumeIdentifier &volIId,
                 std::vector<double> &dims, HepPoint3D &xT) const = 0;

    virtual StatusCode getCorners(const std::vector<double> &dim,
                                  const HepPoint3D &center, 
                                  HepPoint3D *corner) = 0; 

    virtual const AcdUtil::AcdDetectorList& getDetectorList() const = 0;


    virtual bool findDetector(const idents::VolumeIdentifier &volId) const = 0;
};

#endif
