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
#include "AcdUtil/AcdFrameUtil.h"

class AcdGeomMap;
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

static const InterfaceID IID_IAcdGeometrySvc("IAcdGeometrySvc", 1, 4); 

class IAcdGeometrySvc : public virtual IInterface
{
public:

    //! Constructor of this form must be provided

    // Interface ID
    static const InterfaceID& interfaceID() { return IID_IAcdGeometrySvc; }

    
    // Retrieve stored information

    // Simple constants
    virtual int    numTiles()     const = 0;
    virtual int    numRibbons()   const = 0;

    // Lists and Maps of detector elements

    // All the elements
    virtual const AcdUtil::AcdDetectorList& getDetectorList() const = 0;

    // How many sub-elements for a given AcdId
    virtual const std::map<idents::AcdId, int>& getAcdIdVolCountCol() const = 0;

    // Does a volume ID correspond to a known detector elemnt
    virtual bool findDetector(const idents::VolumeIdentifier &volId) const = 0;
    
    // A map from AcdID to the represenations used in AcdRecon
    virtual AcdGeomMap& geomMap() = 0;

    // Get a particular corner Ray
    virtual const Ray getCornerGapRay(unsigned int i) const = 0;


    /// Given an AcdId, provide three vectors of Rays.  Each vector pertains to one set of ribbon segments
    virtual bool fillRibbonData(const idents::AcdId& id,
				std::vector<Ray>& minusSideRays,
				std::vector<Ray>& topRays,
				std::vector<Ray>& plusSideRays, 
				HepTransform3D& minusSideTransform,
				HepTransform3D& topTransform,
				HepTransform3D& plusTransform) = 0;

    // Return half ribbon width
    virtual double ribbonHalfWidth() const = 0;    

    /// Given an AcdId, provide the tile size, center and corners
    virtual bool fillTileData(const idents::AcdId& id, int iVol,
			      HepTransform3D& transform,
			      std::vector<double>& dim, 
			      HepPoint3D& center,
			      HepPoint3D* corner) = 0;

    /// Given an AcdId, provide positions of screw holes in local frame
    virtual bool fillScrewHoleData(const idents::AcdId& id, std::vector< HepPoint3D >& screwHoles) = 0;

    /// Given an AcdId, provide information about which volume edges are shared
    virtual bool fillTileSharedEdgeData(const idents::AcdId& id, 
					const std::vector<double>& dim1, const std::vector<double>& dim2,
					int& sharedEdge1, int& sharedEdge2,
					float& sharedWidth1, float& sharedWidth2) = 0;

    virtual AcdFrameUtil::AcdReferenceFrame getReferenceFrame(const idents::VolumeIdentifier &volId) const = 0;

    virtual StatusCode findCornerGaps() = 0;


};

#endif
