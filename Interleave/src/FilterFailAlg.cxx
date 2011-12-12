/** @file FilterFailAlg.cxx

@brief declaration and definition of the class FilterFailAlg

$Header$

*/

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

/** @class FilterFailAlg
    @brief Gaudi Alg that exists solely to setFilterPassed to false for ending branches
           as part of the interleave algorithm
 */


class FilterFailAlg : public Algorithm {
public:
    FilterFailAlg(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
};

//------------------------------------------------------------------------

//static const AlgFactory<FilterFailAlg>  Factory;
//const IAlgFactory& FilterFailAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(FilterFailAlg);

//------------------------------------------------------------------------
//! ctor
FilterFailAlg::FilterFailAlg(const std::string& name, ISvcLocator* pSvcLocator) : Algorithm(name, pSvcLocator)
{

}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode FilterFailAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode FilterFailAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;

    // We have read the event in, want to now proceed down the interleave branch
    setFilterPassed(false); // since this is on a branch, and we want the sequence to fail

    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FilterFailAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    return sc;
}
