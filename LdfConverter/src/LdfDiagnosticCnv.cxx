// File and Version Information:
//      $Header$
//
// Description:
//      LdfDiagnosticCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfDiagnosticCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfDiagnosticCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"

#include "LdfEvent/DiagnosticData.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfDiagnosticCnv> s_factory;
const ICnvFactory& LdfDiagnosticCnvFactory = s_factory;

LdfDiagnosticCnv::LdfDiagnosticCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/Diagnostic", objType(), "PASS");
}


StatusCode LdfDiagnosticCnv::createObj(IOpaqueAddress* , 
                                    DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF diagnostic data on the TDS
    LdfEvent::DiagnosticData *diag = new LdfEvent::DiagnosticData();
    refpObject = diag;
    
    // Retrieve the LAT data for this event and print out its contents
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    unsigned int i;
    for (i = 0; i < 16; i++) {
        if(myLatData->getTower(i)) {
            if (!myLatData->getTower(i)->getTem().exist()) continue;

            const ldfReader::DiagnosticData* diagnostic = myLatData->getTower(i)->getTem().getDiagnostic();

            if (!diagnostic->exist()) continue;

	    // Added 03.18.2005 by awb to get all Tkr and Cal diagnostics contributions. 
            // Had to hardcode '8' - this should be fixed later.  
	    unsigned int ind;
	    for (ind = 0; ind < 8; ind++) {
	      ldfReader::CalDiagnosticData ldfCalDiag = diagnostic->getCalDiagnosticByIndex(ind);
	      //LdfEvent::CalDiagnosticData cal(diagnostic->getCalDiagnosticByIndex(ind).dataWord());
	      LdfEvent::CalDiagnosticData cal(ldfCalDiag.dataWord(), ldfCalDiag.tower(), ldfCalDiag.layer());
	      diag->addCalDiagnostic(cal);
	    }
	    for (ind = 0; ind < 8; ind++) {
	      ldfReader::TkrDiagnosticData ldfTkrDiag = diagnostic->getTkrDiagnosticByIndex(ind);
	      LdfEvent::TkrDiagnosticData tkr(ldfTkrDiag.dataWord(), ldfTkrDiag.tower(), ldfTkrDiag.gtcc());
	      diag->addTkrDiagnostic(tkr);
	    }

	    // Taken out 03.18.2005 by awb:
            //ldfReader::CalDiagnosticData ldfCalDiag = diagnostic->getCalDiagnostic();
            //LdfEvent::CalDiagnosticData cal(ldfCalDiag.dataWord(), ldfCalDiag.tower(), ldfCalDiag.layer());
            //diag->addCalDiagnostic(cal);

            //ldfReader::TkrDiagnosticData ldfTkrDiag = diagnostic->getTkrDiagnostic();
            // LdfEvent::TkrDiagnosticData tkr(ldfTkrDiag.dataWord(), ldfTkrDiag.tower(), ldfTkrDiag.gtcc());
            // diag->addTkrDiagnostic(tkr);
        }
    }
    return StatusCode::SUCCESS;
}

StatusCode LdfDiagnosticCnv::updateObj(int* , LdfEvent::DiagnosticData* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfDiagnosticCnv");
    log << MSG::DEBUG << "LdfDiagnosticCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}

