// $Id$
// 
//!  \author: Sawyer Gillespie  hgillesp@u.washington.edu
//

#ifndef _H_IGlastDetSvc_
#define _H_IGlastDetSvc_

// includes
#include "Gaudi/Interfaces/IInterface.h"

//forward declarations

class DetectorConverter;
namespace xml {class IFile; }

//!  Access to the Glast detector geometry and IRF I/O

class   IGlastDetSvc : virtual public IInterface {
public:
    
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
};

// IID declaration
extern const IID&   IID_IGlastDetSvc;

#endif  // _H_IGlastDetSvc_
