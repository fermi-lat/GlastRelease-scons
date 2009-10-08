#ifndef AcdFailureModeSvc_H

#define AcdFailureModeSvc_H 1

// $Header$
/** @file
@author H.Kelly
*/


#include "AcdUtil/IAcdFailureModeSvc.h"
#include "GaudiKernel/Service.h"
#include <vector>

/** @class AcdFailureModeSvc
* @brief Service to store and compare to a list of desired failure modes in the ACD.
*
* Author:  H.Kelly
*
*/

class AcdFailureModeSvc : public Service, virtual public IAcdFailureModeSvc {

public:

    AcdFailureModeSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();

    StatusCode finalize();


    /// look for AcdId in list of dead tiles

    bool matchAcdId(idents::AcdId id);


    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);


    static const InterfaceID& interfaceID() {
        return IAcdFailureModeSvc::interfaceID(); 
    }


    /// return the service type
    const InterfaceID&  type () const {
        return IID_IAcdFailureModeSvc;
    }


protected:

    /// process the input list of tiles
    void processDetectorList();


private:

    /// List of detectors from jobOptions
    UnsignedIntegerArrayProperty m_detectorListProperty;

    /// vector of detectors to fail
    std::vector<unsigned int> m_detectorList;


};

#endif // AcdFailureModeSvc_H

