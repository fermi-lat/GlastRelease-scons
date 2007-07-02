
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
#include <map>

class AcdGeometrySvc : public Service,
        virtual public IAcdGeometrySvc
{
public:
    
    AcdGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~AcdGeometrySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return IAcdGeometrySvc::interfaceID(); 
    }
    /// return the service type
    const InterfaceID& type() const;

    StatusCode getConstants();
    StatusCode getDetectorListFromGeometry();

    const AcdUtil::AcdDetectorList& getDetectorList() const { 
        return m_detectorCol; };

    const std::map<idents::AcdId, int>& getAcdIdVolCountCol() const { 
        return m_acdId_volCount; };

    bool findDetector(const idents::VolumeIdentifier &volId) const;

    int numTiles() const { return m_numTiles; };
    int numRibbons() const {return m_numRibbons; };
  
    StatusCode getCorners(const std::vector<double> &dim,
                          const HepPoint3D &center, HepPoint3D *corner);

    StatusCode getDimensions(const idents::VolumeIdentifier &volIId, 
                                    std::vector<double> &dims, 
                                    HepPoint3D &xT) const;

    StatusCode getDetectorDimensions(const idents::VolumeIdentifier &volIId, 
                                    std::vector<double> &dims, 
                                    HepPoint3D &xT) const;

    StatusCode findCornerGaps();
    const Ray getCornerGapRay(unsigned int index) const;

    /// Given an AcdId, provide three vectors of Rays.  Each vector pertains to one set of ribbon segments
    bool fillRibbonRays(const idents::AcdId& id,
                 std::vector<Ray>& minusSideRays,
                 std::vector<Ray>& topRays,
                 std::vector<Ray>& plusSideRays, bool increasing = true);

    /// Given an AcdId for a ribbon, provide the transformation to the center of each set of ribbon segments
    virtual bool fillRibbonTransform(int face,
				     const Ray& ribbon,
				     HepTransform3D& transform);

    virtual double ribbonHalfWidth() const;

    /// Given an AcdId, provide the tile size, center and corners
    virtual bool fillTileData(const idents::AcdId& id, int iVol,
			      std::vector<double>& dim, 
			      HepPoint3D& center,
			      HepPoint3D* corner);

    /// Given an AcdId, provide transform to tile frame
    virtual bool fillTileTransform(const idents::AcdId& id, int iVol,
				   HepTransform3D& transform);

    /// Given an AcdId, provide positions of screw holes in local frame
    virtual bool fillScrewHoleData(const idents::AcdId& id, std::vector< HepPoint3D >& screwHoles);

    /// Given an AcdId, provide information about which volume edges are shared
    virtual bool fillTileSharedEdgeData(const idents::AcdId& id, 
					const std::vector<double>& dim1, const std::vector<double>& dim1,
					int& sharedEdge1, int& sharedEdge2,
					float& sharedWidth1, float& sharedWidth2);


private:

    void clear();
    
    /// pointer to the detector service
    IGlastDetSvc *m_glastDetSvc;

    int m_numTiles, m_numRibbons;
    int m_numXtowers, m_numYtowers;

    Point m_cornerGapStartPoint[4], m_cornerGapEndPoint[4];
    Vector m_cornerGapVec[4];

    AcdUtil::AcdDetectorList m_detectorCol;

    /// A count of volumes associated with each AcdId
    std::map<idents::AcdId, int> m_acdId_volCount;
};

#endif // ACDGEOMETRYSVC_H
