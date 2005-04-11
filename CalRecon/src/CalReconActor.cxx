
#include "CalReconActor.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IToolSvc.h"

CalReconActor::CalReconActor()
 : m_kernel(0)
 {}

CalReconActor::~CalReconActor()
 {}

StatusCode CalReconActor::initialize( ISvcLocator * svcLocator )
 {
  IMessageSvc * msgSvc ;
  svcLocator->service("MessageSvc",msgSvc,true) ;
  MsgStream log(msgSvc, "CalReconActor") ;
    
  IToolSvc * toolSvc ;
  if ((svcLocator->service("ToolSvc",toolSvc,true)).isFailure())
   {
    log<<MSG::ERROR<<"Unable to find the ToolSvc"<<endreq ;
    return StatusCode::FAILURE ;
   }
   
  if ((toolSvc->retrieveTool("CalReconKernel",m_kernel)).isFailure())
   {
    log<<MSG::ERROR<<"Unable to retrieve the kernel"<<endreq ;
    return StatusCode::FAILURE ;
   }
      
  return StatusCode::SUCCESS ; 
 }
    
CalReconKernel * CalReconActor::getKernel()
 { return m_kernel ; }

