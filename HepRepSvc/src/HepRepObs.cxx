/** @file HepRepObs.cxx 
   @brief Implementation file for HepRepObs 


   finds instances of IRegister tools

 $Header$

*/

#include "HepRepObs.h"
#include "HepRepSvc/IRegister.h"

#include <iterator>
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataStoreItem.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/SmartDataPtr.h"


HepRepObs::HepRepObs():IToolSvc::Observer()
{
   
}


void HepRepObs::onCreate(IAlgTool& tool) {


    IRegister* gtool;
    StatusCode status =tool.queryInterface( IRegister::interfaceID(), (void**)&gtool);
    if( status.isSuccess() ){
        gtool->registerMe(m_hepRepSvc);
    }


}


