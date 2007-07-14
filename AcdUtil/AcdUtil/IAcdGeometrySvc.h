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
#include <vector>

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

    /// Given an AcdId, provide three vectors of Rays.  Each vector pertains to one set of ribbon segments
    virtual bool fillRibbonRays(const idents::AcdId& id,
                 std::vector<Ray>& minusSideRays,
                 std::vector<Ray>& topRays,
                 std::vector<Ray>& plusSideRays, bool increasing = true) = 0;

    /// Given an AcdId for a ribbon, provide the transformation to the center of each set of ribbon segments
    virtual bool fillRibbonTransforms(const idents::AcdId& id,
				      HepTransform3D& minusSideTransform,
				      HepTransform3D& topTransform,
				      HepTransform3D& plusTransform) = 0;		    

    virtual double ribbonHalfWidth() const = 0;

    /// Given an AcdId, provide the tile size, center and corners
    virtual bool fillTileData(const idents::AcdId& id, int iVol,
			      int& face,
			      std::vector<double>& dim, 
			      HepPoint3D& center,
			      HepPoint3D* corner) = 0;

    /// Given an AcdId, provide transform to tile frame
    virtual bool fillTileTransform(const idents::AcdId& id, int iVol,
				   HepTransform3D& transform) = 0;

    /// Given an AcdId, provide positions of screw holes in local frame
    virtual bool fillScrewHoleData(const idents::AcdId& id, std::vector< HepPoint3D >& screwHoles) = 0;

    /// Given an AcdId, provide information about which volume edges are shared
    virtual bool fillTileSharedEdgeData(const idents::AcdId& id, 
					const std::vector<double>& dim1, const std::vector<double>& dim1,
					int& sharedEdge1, int& sharedEdge2,
					float& sharedWidth1, float& sharedWidth2) = 0;

};

#endif
