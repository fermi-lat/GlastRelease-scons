// File and Version Information:
//      $Header$
//
// Description:
//      LdfCalDigiCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfCalDigiCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfCalDigiCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"

#include "Event/Digi/CalDigi.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfCalDigiCnv> s_factory;
const ICnvFactory& LdfCalDigiCnvFactory = s_factory;

LdfCalDigiCnv::LdfCalDigiCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/Digi/CalDigiCol", objType(), "PASS");
}


StatusCode LdfCalDigiCnv::createObj(IOpaqueAddress* , DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    //refpObject = new Event::EventHeader();
    Event::CalDigiCol *digiCol = new Event::CalDigiCol;
    refpObject = digiCol;

    
    // Retrieve the LAT data for this event and print out its contents
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    unsigned int iTower;
    for (iTower = 0; iTower < 16; iTower++) {
       ldfReader::TowerData *tem = myLatData->getTower(iTower);
        if (!tem) {
            //printf("Tem does not exist\n");
        } else {
            const std::map<unsigned int, ldfReader::CalDigi*> calData = tem->getCalDigiCol();
            std::map<unsigned int, ldfReader::CalDigi*>::const_iterator it;
            for (it = calData.begin(); it != calData.end(); it++) {
                int layer = it->second->getLayerNum();
                int col = it->second->getColumnNum();
                idents::CalXtalId xtalId(iTower, layer, col);
                ldfReader::CalDigi::CalTrigMode mode = it->second->getMode();
                if (mode == ldfReader::CalDigi::BESTRANGE) {
                    const ldfReader::CalDigi::CalReadout *readout = it->second->getXtalReadout(0);
                    char rangeP = readout->getRange(ldfReader::CalDigi::POS);
                    unsigned short adcP = readout->getAdc(ldfReader::CalDigi::POS);
                    char rangeM = readout->getRange(ldfReader::CalDigi::NEG);
                    unsigned short adcM = readout->getAdc(ldfReader::CalDigi::NEG);
                    Event::CalDigi::CalXtalReadout read = Event::CalDigi::CalXtalReadout(rangeP, adcP, rangeM, adcM);
                    Event::CalDigi* curDigi = new Event::CalDigi(idents::CalXtalId::BESTRANGE, xtalId);
                    curDigi->addReadout(read);
                    digiCol->push_back(curDigi);

                } else { // ALLRANGE
                    Event::CalDigi* curDigi = new Event::CalDigi(idents::CalXtalId::ALLRANGE, xtalId);
                    unsigned int j;
                    for (j=0; j<4; j++) {
                        const ldfReader::CalDigi::CalReadout *readout = it->second->getXtalReadout(j);
                        char rangeP = readout->getRange(ldfReader::CalDigi::POS);
                        unsigned short adcP = readout->getAdc(ldfReader::CalDigi::POS);
                        char rangeM = readout->getRange(ldfReader::CalDigi::NEG);
                        unsigned short adcM = readout->getAdc(ldfReader::CalDigi::NEG);
                        Event::CalDigi::CalXtalReadout read = Event::CalDigi::CalXtalReadout(rangeP, adcP, rangeM, adcM);
                        curDigi->addReadout(read);
                    }
                    digiCol->push_back(curDigi);
                }
            }

        }
    }
    return StatusCode::SUCCESS;
}

StatusCode LdfCalDigiCnv::updateObj(int* , Event::EventHeader* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfCalDigiCnv");
    log << MSG::DEBUG << "LdfCalDigiCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}

