// $Id$
// 
//!  \author: Sawyer Gillespie  hgillesp@u.washington.edu
//

#ifndef _H_IGlastDetSvc_
#define _H_IGlastDetSvc_

// includes
#include "GaudiKernel/IInterface.h"

//forward declarations
class DMmanager;
class DetectorConverter;
class GlastDetector;
class IGeometry;
class IMedia;
class HepTransform3D;
namespace xml {class IFile; }
namespace idents{class VolumeIdentifier;}

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IGlastDetSvc(901, 1 , 0); 
//!  Access to the Glast detector geometry and IRF I/O

class   IGlastDetSvc : virtual public IInterface {
public:
  
    //! detModel interface to retrive numeric constants
    virtual StatusCode getNumericConstByName(std::string, double*)=0;

    //! new detModel interface, will call back. 
    virtual void accept(IGeometry& geom)=0;

    //! new detModel interface, will call back. 
    virtual void accept(IMedia& media)=0;

    /// retrive the 3D transformation of a volume given a valid ID
    virtual StatusCode getTransform3DByID(idents::VolumeIdentifier,HepTransform3D*)=0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IGlastDetSvc; }

};



#endif  // _H_IGlastDetSvc_
