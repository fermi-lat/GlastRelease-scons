
#include "PropertiesCheckAlg.h"
#include <GaudiKernel/AlgFactory.h>
#include <GaudiKernel/IService.h>
#include <GaudiKernel/IJobOptionsSvc.h>
#include <GaudiKernel/DeclareFactoryEntries.h>
#include <GaudiKernel/MsgStream.h>
#include <vector>

DECLARE_ALGORITHM_FACTORY(PropertiesCheckAlg) ;

PropertiesCheckAlg::PropertiesCheckAlg( const std::string & name, ISvcLocator * locator )
 : Algorithm(name,locator)
 {
  declareProperty("targets", m_targets) ;
 }

StatusCode PropertiesCheckAlg::initialize()
 {
  MsgStream log(msgSvc(),name()) ;
  IService * joSvc =0 ;
  IJobOptionsSvc * joInt =0 ;
  StatusCode joSc ;
  joSc = service("JobOptionsSvc",joSvc,true) ;
  if (joSc.isSuccess())
   { joSc = joSvc->queryInterface(IID_IJobOptionsSvc,(void**)&joInt) ; }
  if (joSc.isSuccess())
   {
    std::vector< std::string > clients = joInt->getClients() ;
    std::vector< std::string >::iterator client ;
    for ( client = clients.begin() ; client != clients.end() ; ++client )
     {
      const std::vector< std::string > & targets = m_targets ;
      std::vector< std::string >::iterator target ;
      for ( target = targets.begin() ; target != targets.end() ; ++target )
       {
        if ((*target)==(*client))
         {
          const std::vector< const Property * > * properties = joInt->getProperties(*client) ;  
          std::vector< const Property * >::const_iterator property ;
          for ( property = properties->begin() ; property != properties->end() ; ++property )
           {
            Property * p = const_cast<Property *>(*property) ;
            p->declareReadHandler(&PropertiesCheckAlg::readPropertyCallBack,this) ;
            m_properties[p->name()] = 0 ;
           }
          break ;
         }
       }
     }
   }  
  return StatusCode::SUCCESS ;
 }

void PropertiesCheckAlg::readPropertyCallBack( Property & p )
 { m_properties[p.name()]++ ; }

StatusCode PropertiesCheckAlg::finalize()
 {
  MsgStream log(msgSvc(),name()) ;
  StatusCode sc = StatusCode::SUCCESS ;
  std::map< std::string, int >::iterator p ;
  for ( p = m_properties.begin() ; p != m_properties.end() ; ++p )
   {
    if (((*p).second)==0)
     {
      log<<MSG::ERROR<<(*p).first<<" NOT USED"<<endreq ;
      sc = StatusCode::FAILURE ;
     }
   }
  return sc ;
 }


