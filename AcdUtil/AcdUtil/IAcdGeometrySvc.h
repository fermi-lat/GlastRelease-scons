/** @file IAcdGeometrySvc.h
 @brief Abstract interface to AcdGeometrySvc (q.v.)

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
class AcdRibbonSegment;

#include <string>
#include <map>
#include <vector>

/** 
 * @class IAcdGeometrySvc
 *
 * @brief Interface to ACD geometry.
 *
 * Most clients just use the AcdGeomMap to get the descriptions of individual elements.
 * However there are also some other useful functions:
 * - Used in filling the AcdGeomMap
 *   - getting the number of elements, and a list of all the elements
 *   - getting the number of volumes that make up a particular elements
 *   - functions to fill the AcdTileDim, and AcdRibbonDim data members
 *   - determining the intermediate reference frame of a particular volume of an ACD element
 * - Used in reconstruction
 *   - getting the rays associated with the corner gaps of the ACD
 *   - deciding if an idents::VolumeIdentifier is an ACD detector element or not
 * 
 * @author Heather Kelly 
 */

static const InterfaceID IID_IAcdGeometrySvc("IAcdGeometrySvc", 1, 5); 

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

    /// All the elements
    virtual const AcdUtil::AcdDetectorList& getDetectorList() const = 0;

    /// How many sub-elements for a given AcdId
    virtual const std::map<idents::AcdId, int>& getAcdIdVolCountCol() const = 0;

    /// Does a volume ID correspond to a known detector elemnt
    virtual bool findDetector(const idents::VolumeIdentifier &volId) const = 0;
    
    /// A map from AcdID to the represenations used in AcdRecon
    virtual AcdGeomMap& geomMap() = 0;

    /// Get a particular corner Ray
    virtual const Ray getCornerGapRay(unsigned int i) const = 0;

    /// Given an AcdId, provide three vectors of Rays.  Each vector pertains to one set of ribbon segments
    virtual bool fillRibbonData(const idents::AcdId& id,
				std::vector<AcdRibbonSegment*>& segs,
				int& topIdx, int& plusIdx) = 0;

    // Return half ribbon width
    virtual double ribbonHalfWidth() const = 0;    

    /// Given an AcdId, provide the tile size, center and corners
    virtual bool fillTileData(const idents::AcdId& id, int iVol,
			      HepGeom::Transform3D& transform,
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

    /// Given a VolumeIdentifier, provide the intermediate reference frame
    virtual AcdFrameUtil::AcdReferenceFrame getReferenceFrame(const idents::VolumeIdentifier &volId) const = 0;

    virtual StatusCode findCornerGaps() = 0;

    virtual StatusCode getNextTileCorners(const idents::AcdId& id, int dir, 
					  HepPoint3D& c1, HepPoint3D& c2, bool& isRealGap) = 0;


};

#endif
