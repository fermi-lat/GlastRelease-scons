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
namespace xml {class IFile; }

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
#if 0
    //! open an IRF file
    virtual StatusCode openIRF(std::string filename)=0;
    
    //! read the next event from the IRF file return fail if EOF
    virtual StatusCode  readIRF ()=0;

    
    //! return the number of detectors currently with data
    virtual int detectors_with_data()const=0;
  
    //! accept a visitor to traverse the structure note that must be on a const object
    virtual void accept(DetectorConverter&)const=0;

    //! access to the IFile containing detector constants
    virtual const xml::IFile* iniFile()const=0;

    //! set new root detector
    virtual void setDetector(GlastDetector* d)=0;
#endif
	    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IGlastDetSvc; }

};



#endif  // _H_IGlastDetSvc_
