
#ifndef ACDGEOMETRYSVC_H
#define ACDGEOMETRYSVC_H 1

/** 
 * @class AcdGeometrySvc
 *
 * @brief Supplies the geometry constants and calculations for ACD 
 *
 * The constants all flow from GlastDetSvc.
 * 
 * @author Heather Kelly 
 *
 * $Header$
 */

#include "GaudiKernel/Service.h"

#include "AcdUtil/IAcdGeometrySvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"
#include "idents/AcdId.h"

#include "geometry/Ray.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "AcdUtil/AcdDetectorList.h"
#include "AcdUtil/AcdGeomMap.h"
#include <map>

class AcdGeometrySvc : public Service,
        virtual public IAcdGeometrySvc
{
public:
    
    AcdGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~AcdGeometrySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    

    /// return the Interface ID
    static const InterfaceID& interfaceID() {
        return IAcdGeometrySvc::interfaceID(); 
    }
    /// return the service type
    const InterfaceID& type() const;


    /// Return simple data

    /// Simple constants
    int numTiles() const { return m_numTiles; };
    int numRibbons() const {return m_numRibbons; };
   
    /// query by name
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

    /// get the list of detector elements 
    const AcdUtil::AcdDetectorList& getDetectorList() const { 
        return m_detectorCol; };

    /// How many sub-elements for a given AcdId
    const std::map<idents::AcdId, int>& getAcdIdVolCountCol() const { 
        return m_acdId_volCount; };

    /// Does a volume ID correspond to a known detector elemnt
    bool findDetector(const idents::VolumeIdentifier &volId) const;

    /// map from AcdID to the represenations used in AcdRecon
    AcdGeomMap& geomMap() {
        return m_geomMap; }

    /// Get a particular corner Ray
    const Ray getCornerGapRay(unsigned int index) const;

    /// Given an AcdId, provide three vectors of Rays.  Each vector pertains to one set of ribbon segments
    bool fillRibbonRays(const idents::AcdId& id,
			std::vector<Ray>& minusSideRays,
			std::vector<Ray>& topRays,
			std::vector<Ray>& plusSideRays, bool increasing = true);
    
    /// Given an AcdId for a ribbon, provide the transformation to the center of each set of ribbon segments
    virtual bool fillRibbonTransforms(const idents::AcdId& id,
				      HepTransform3D& minusSideTransform,
				      HepTransform3D& topTransform,
				      HepTransform3D& plusTransform);	

    /// Return half ribbon width
    virtual double ribbonHalfWidth() const;

 
    /// Given an AcdId, provide the tile size, center and corners
    virtual bool fillTileData(const idents::AcdId& id, int iVol,
			      HepTransform3D& transform,
			      std::vector<double>& dim, 
			      HepPoint3D& center,
			      HepPoint3D* corner);

    /// Given an AcdId, provide positions of screw holes in local frame
    virtual bool fillScrewHoleData(const idents::AcdId& id, std::vector< HepPoint3D >& screwHoles);

    /// Given an AcdId, provide information about which volume edges are shared
    virtual bool fillTileSharedEdgeData(const idents::AcdId& id, 
					const std::vector<double>& dim1, const std::vector<double>& dim2,
					int& sharedEdge1, int& sharedEdge2,
					float& sharedWidth1, float& sharedWidth2);

    AcdFrameUtil::AcdReferenceFrame getReferenceFrame(const idents::VolumeIdentifier &volId);

protected:

    /// Used in initializing 
    StatusCode getConstants();
    StatusCode getDetectorListFromGeometry();
    StatusCode findCornerGaps();

    // Utilities
    StatusCode getDimensions(const idents::VolumeIdentifier &volIId, 
			     std::vector<double> &dims, 
			     HepPoint3D &xT) const;
  
    StatusCode getDetectorDimensions(const idents::VolumeIdentifier &volIId, 
				     std::vector<double> &dims, 
				     HepPoint3D &xT) const;
  

private:

    void clear();

    typedef std::vector<idents::VolumeIdentifier>::const_iterator VolIdIter;
    
    /// pointer to the detector service
    IGlastDetSvc *m_glastDetSvc;

    AcdGeomMap m_geomMap;

    int m_numTiles, m_numRibbons;
    int m_numXtowers, m_numYtowers;

    unsigned m_eLATACD;
    unsigned m_eACDTile,  m_eACDRibbon;
    unsigned m_eMeasureX,  m_eMeasureY;
    unsigned m_eACDTopFace, m_eACDXNegFace, m_eACDYNegFace;
    unsigned m_eACDXPosFace, m_eACDYPosFace;

    double m_ribbonHalfWidth;

    Point m_cornerGapStartPoint[4], m_cornerGapEndPoint[4];
    Vector m_cornerGapVec[4];

    AcdUtil::AcdDetectorList m_detectorCol;

    /// A count of volumes associated with each AcdId
    std::map<idents::AcdId, int> m_acdId_volCount;
};

#endif // ACDGEOMETRYSVC_H
