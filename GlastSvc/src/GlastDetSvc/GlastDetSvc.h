// $Header$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastDetSvc_
#define _H_GlastDetSvc_

// includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IEvtSelector.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// forward declarations
class Instrument;
class GlastDetector;
namespace xml { class IFile; }
class DMmanager;
class HepTransform3D;
class HepPoint3D;
namespace idents{class VolumeIdentifier;}

/*!  This Gaudi service provides access to a tree of GlastDetector objects.
*/
class GlastDetSvc : public Service, 
    virtual public IGlastDetSvc
{
public:
    
    //! Create the service
    //! @param name The name of the service
    //!
    GlastDetSvc ( const std::string& name, ISvcLocator* al );
    
    /// destructor
    virtual ~GlastDetSvc ();

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);
    
    /// return the service type
    const IID& type() const;
    
    /// perform initializations for this service. This will
    ///    read in and create the GlastDetector hierarchy from the file
    ///    specified by the PersistencyFile property, unless it is
    ///    '(none)'
    StatusCode initialize ();
    
    /// kill off the detector hierarchy if it exists
    StatusCode finalize ();
    //! new detModel interface, will call back. 
    virtual void accept(IMedia& media);

    //! start a visitor of the detModel geometry description (implements IGlastDetSvc)
    virtual void accept(IGeometry& geom);

    //! detModel interface to retrieve numeric constants
    virtual StatusCode getNumericConstByName(std::string, double*);

    /// retrive the 3D transformation of a volume given a valid ID
    virtual StatusCode getTransform3DByID(idents::VolumeIdentifier,HepTransform3D*);

    /// retrive the type and dimensions of a volume given a valid ID
    virtual StatusCode  getShapeByID(idents::VolumeIdentifier id,
                                     std::string*, 
                                     std::vector<double>*);

    
    /// compute strip id from local coordinate for a tkr plane
    virtual unsigned int stripId (double x);
    virtual unsigned int stripId (const HepPoint3D&  x);
    /// (-) if non-active or gap (P is in local system)    
    virtual double insideActiveArea (const HepPoint3D& p);  

    /// location of strip ix in local coords
    virtual double stripLocalX ( unsigned int istrip); 
	
    /// location of strip ix in local coords
    virtual double stripLocalXDouble ( double strip);   
	
	/// return plane coordinate given wafer coordinate
	virtual HepPoint3D siPlaneCoord( const HepPoint3D &p, idents::VolumeIdentifier id);


private:
    
    // Data Members
    Instrument*		m_instrument;
    std::string		m_xmlFile;
    std::string		m_iniFile;
    DMmanager*          m_dm;
    std::string         m_xmlfile;
    std::string         m_topvol;
    std::string         m_visitorMode;
};

#endif // _H_GlastDetSvc
