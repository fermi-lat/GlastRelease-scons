
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

#include "AcdUtil/AcdDetectorList.h"

class AcdGeometrySvc : public Service,
        virtual public IAcdGeometrySvc
{
public:
    
    AcdGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~AcdGeometrySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return IAcdGeometrySvc::interfaceID(); 
    }
    /// return the service type
    const IID& type() const;

    StatusCode getConstants();
    StatusCode getDetectorListFromGeometry();

    const AcdUtil::AcdDetectorList& getDetectorList() const { 
        return m_detectorCol; };

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

private:

    void clear();
    
    /// pointer to the detector service
    IGlastDetSvc *m_glastDetSvc;

    int m_numTiles, m_numRibbons;
    int m_numXtowers, m_numYtowers;

    AcdUtil::AcdDetectorList m_detectorCol;
};

#endif // ACDGEOMETRYSVC_H
