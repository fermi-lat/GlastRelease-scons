/** @file IAcdGeometrySvc.h
 @brief Abstract interface to TkrGeometrySvc (q.v.)

  $Header$
*/

#ifndef __IACDGEOMETRYSVC_H
#define __IACDGEOMETRYSVC_H 1

#include "GaudiKernel/IInterface.h"


#include "CLHEP/Geometry/Point3D.h"

// TU: Hacks for CLHEP 1.9.2.2 and beyond
#ifndef HepPoint3D
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "geometry/Ray.h"

#include "AcdUtil/AcdDetectorList.h"

#include <string>
#include <map>

/** 
 * @class IAcdGeometrySvc
 *
 * @brief Abstract interface to AcdGeometrySvc 
 * 
 * @author Heather Kelly 
 */

static const InterfaceID IID_IAcdGeometrySvc("IAcdGeometrySvc", 1, 2); 

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

    virtual const std::map<idents::AcdId, int>& getAcdIdVolCountCol() const = 0;

    virtual bool findDetector(const idents::VolumeIdentifier &volId) const = 0;

    virtual StatusCode findCornerGaps() = 0;
    virtual const Ray getCornerGapRay(unsigned int i) const = 0;

};

#endif
