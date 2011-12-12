/** @file Count.cxx

*/
#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataSvc.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <map>
/** @class Count
 * @brief Alg that counts and records the value in the jobinfo tuple
 * 
 * $Header$
*/
class Count : public Algorithm
{   
public:
    Count(const std::string& name, ISvcLocator* pSvcLocator); 
    ~Count();
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
private:
    INTupleWriterSvc* m_ntupleWriterSvc;
    // have to save the counts in a static, since this object gets deleted before the tuple is filled :-(
    static std::map<std::string, int> s_map; 
};

//  factory stuff
//static const AlgFactory<Count>  Factory;
//const IAlgFactory& CountFactory = Factory;
DECLARE_ALGORITHM_FACTORY(Count);

std::map<std::string,int> Count::s_map;

Count::Count(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator) 
{
    s_map[name]=0;
}
Count::~Count(){
}

StatusCode Count::initialize() { 
    MsgStream log(msgSvc(), name());
    StatusCode sc =StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    sc = service("RootTupleSvc", m_ntupleWriterSvc);
    
    if(sc.isFailure())    {
        log << MSG::ERROR << "Could not locate the ntupleSvc" <<endreq;
    }
    return sc;   
}


StatusCode Count::execute() {

    ++ s_map[name()];
    return StatusCode::SUCCESS;
}


StatusCode Count::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "counted " << s_map[name()] << " times" << endreq;
    std::map<std::string, int>::iterator p = s_map.find(name());
    int * i = &(p->second);
    m_ntupleWriterSvc->addItem("jobinfo", name(), i);
    return StatusCode::SUCCESS;
}

