
#ifndef _H_IGlastDetSvc_
#define _H_IGlastDetSvc_

#include "GaudiKernel/IInterface.h"
#include <vector>
#include <string>

class DMmanager;
class DetectorConverter;
class GlastDetector;
class IGeometry;
class IMedia;
class HepTransform3D;
class HepPoint3D;
class IDmapBuilder;

namespace idents{class VolumeIdentifier;}

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IGlastDetSvc(901, 1 , 0); 
//!  Access to the Glast detector geometry and IRF I/O

/** @class IGlastDetSvc
* @brief 
*
* @author Sawyer Gillespie
*
* $Header$
*/

class   IGlastDetSvc : virtual public IInterface {
public:
  
    /// detModel interface to retrive numeric constants
    virtual StatusCode getNumericConstByName(std::string, double*)=0;

    /** detModel interface to retrieve integer numeric constants.
        Returns FAILURE and stores no value in val if constant doesn't
        exist or has not been declared of type int in the xml source.
    */
    virtual StatusCode getNumericConstByName(std::string, int*)=0;

    /// Return Volume identifer of top volume relative to world
    virtual idents::VolumeIdentifier getIDPrefix()=0;


    //! new detModel interface, will call back. 
    virtual void accept(IGeometry& geom)=0;

    //! new detModel interface, will call back. 
    virtual void accept(IMedia& media)=0;

    /// retrieve the 3D transformation of a volume given a valid ID
    virtual StatusCode getTransform3DByID(idents::VolumeIdentifier,HepTransform3D*)=0;

    /// Return transform of top volume relative to world
    virtual const HepTransform3D& getTransform3DPrefix()=0;    

    /// retrive the type and dimensions of a volume given a valid ID
    virtual StatusCode  getShapeByID(idents::VolumeIdentifier id,
                                     std::string*, 
                                     std::vector<double>*) = 0;


    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IGlastDetSvc; }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//            SiliconPlaneGeometry interface calls
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// compute strip id from local coordinate for a tkr plane
    virtual unsigned int stripId (double x)=0;

    /// location of strip ix in local coords
    virtual double stripLocalX( unsigned int stripid)=0; 
	
    /// location of strip ix in local coords
    virtual double stripLocalXDouble( double stripid)=0;    

	/// returns plane coord given ladder coord
	virtual HepPoint3D siPlaneCoord( const HepPoint3D& p, idents::VolumeIdentifier id)=0;

    /// (-) if non-active or gap (P is in local system)    
    virtual double insideActiveArea (const HepPoint3D& p)=0;  

};



#endif  // _H_IGlastDetSvc_
