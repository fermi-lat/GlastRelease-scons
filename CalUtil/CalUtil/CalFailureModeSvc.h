#ifndef CalFailureModeSvc_H
#define CalFailureModeSvc_H 1

// Include files
#include "ICalFailureModeSvc.h"
#include "GaudiKernel/Service.h"

/** @class CalFailureModeSvc
* @brief Service to store and compare to a list of desired failure modes in
* the CAL.
*
* Author:  R.Dubois
*
*/

class CalFailureModeSvc : public Service, virtual public ICalFailureModeSvc {
    
public:

    CalFailureModeSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    /// get the list of enabled failure mode conditions
    int getFailureConditions() {return m_failureModes;};

    enum {TOWER, TOWERAFEE, TOWERCONTROL};

    /// look for crystal in list of dead towers
    bool matchChannel(idents::CalXtalId id, idents::CalXtalId::XtalFace face);

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ICalFailureModeSvc::interfaceID(); 
    }

    /// return the service type
    const IID& type() const;

protected:

    /// look for crystal in list of dead towers
    bool matchTower(idents::CalXtalId id);

    /// look for crystal in list of dead AFEE boards
    bool matchTowerAfee(idents::CalXtalId id, idents::CalXtalId::XtalFace face);

    /// look for crystal in list of dead controllers
    bool matchTowerController(idents::CalXtalId id, idents::CalXtalId::XtalFace face);

    /// process the input list of towers
    void processTowerList();

    /// process the input list of tower/layer pairs
    void processTowerAfeeList();

    /// process the input list of tower/cable controller pairs
    void processTowerControllerList();

private:

    /// List of towers from jobOptions
    StringArrayProperty m_towerListProperty;

    /// List of towers,layers from jobOptions
    StringArrayProperty m_towerAfeeListProperty;

    /// List of towers,layers from jobOptions
    StringArrayProperty m_towerControllerListProperty;

    /// bitmap of failure modes
    int m_failureModes;

    /// vector of towers to fail
    std::vector<int> m_towerList;

    /// map of towers/layers to fail
    std::map<int,std::vector<int> > m_towerAfeeList;

    /// map of towers/layers to fail
    std::map<int,std::vector<int> > m_towerControllerList;
};

#endif // CalFailureModeSvc_H

