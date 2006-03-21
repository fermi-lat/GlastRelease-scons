/** @file IGlastDetSvc.h
* @brief Interface to the GlastDetSvc (q.v.)
*
* @author Sawyer Gillespie
* @author Leon Rochester
*
* $Header$
*/


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
namespace HepGeom {
    class Transform3D;
    template <class T> class Point3D;
}
// TU: Hack for CLHEP 1.9.2.2 
typedef HepGeom::Point3D<double> HepPoint3D;
typedef HepGeom::Transform3D     HepTransform3D;
class IDmapBuilder;

namespace idents{class VolumeIdentifier;}

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IGlastDetSvc(901, 3 , 0); 
///  Access to the Glast detector geometry

/** @class IGlastDetSvc
* @brief Interface to GlastDetSvc
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

    /// new detModel interface, will call back. 
    virtual void accept(IGeometry& geom)=0;

    /// new detModel interface, will call back. 
    virtual void accept(IMedia& media)=0;

    /// retrieve the 3D transformation of a volume given a valid ID
    virtual StatusCode getTransform3DByID(idents::VolumeIdentifier,HepTransform3D*)=0;

    /// Return transform of top volume relative to world
    virtual const HepGeom::Transform3D& getTransform3DPrefix()=0;    

    /// retrieve the type and dimensions of a volume given a valid ID
    virtual StatusCode  getShapeByID(idents::VolumeIdentifier id,
                                     std::string*, 
                                     std::vector<double>*) = 0;
    /// retrieve name of top volume
    virtual std::string getTopVolume() = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IGlastDetSvc; }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//            SiliconPlaneGeometry interface calls
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// compute strip id from local coordinate for a tkr plane
    virtual unsigned int stripId (double x)=0;

    /// location of strip ix in local coords
    virtual double stripLocalX( double stripid)=0; 
    
    /// returns plane coord given ladder coord
    virtual HepPoint3D siPlaneCoord( const HepPoint3D& p, idents::VolumeIdentifier id)=0;

    /// (-) if non-active or gap (P is in local system)    
    virtual double insideActiveArea (const HepPoint3D& p)=0; 

    // methods to return separate x and y values for insideActive
    /// returns distance to x boundary, same convention as InsideActiveArea()
    virtual double insideActiveLocalX( const HepPoint3D& p)=0;
    /// returns distance to y boundary, same convention as InsideActiveArea()
    virtual double insideActiveLocalY( const HepPoint3D& p)=0;

    
    /// return the global position of a strip in the plane, will accept int or double
    virtual HepPoint3D getStripPosition(idents::VolumeIdentifier volId, double stripid)=0;
};

#endif  // _H_IGlastDetSvc_
