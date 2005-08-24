// File and Version Information:
//      $Header$
//
// Description:
//      LdfAcdDigiCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfAcdDigiCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfAcdDigiCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"

//#include "Event/Digi/AcdDigi.h"
#include "idents/AcdId.h"

#include <map>

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfAcdDigiCnv> s_factory;
const ICnvFactory& LdfAcdDigiCnvFactory = s_factory;

LdfAcdDigiCnv::LdfAcdDigiCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event/Digi/AcdDigiCol path on the TDS.
    declareObject("/Event/Digi/AcdDigiCol", objType(), "PASS");
}


StatusCode LdfAcdDigiCnv::createObj(IOpaqueAddress* , DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    MsgStream log(msgSvc(), "LdfAcdDigiCnv");
    Event::AcdDigiCol *digiCol = new Event::AcdDigiCol;
    refpObject = digiCol;

    // Retrieve the LAT data for this event and grab the ACD digi data
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    if (!myLatData) return StatusCode::SUCCESS;
    const std::map<const char*, ldfReader::AcdDigi*> acdCol = myLatData->getAcdCol();
    std::map<const char*, ldfReader::AcdDigi*>::const_iterator thisAcdDigi;
    
    for (thisAcdDigi = acdCol.begin(); thisAcdDigi != acdCol.end(); thisAcdDigi++) {
        const char *tileName = thisAcdDigi->second->getTileName();
        int tileNumber = thisAcdDigi->second->getTileNumber();
        unsigned int tileId = thisAcdDigi->second->getTileId();
		short layer, face, row, col;
		base10ToAcdId(tileId, layer, face, row, col);
        idents::AcdId identsId(layer, face, row, col);

        const std::vector<ldfReader::AcdDigi::AcdPmt> readoutCol = thisAcdDigi->second->getReadout();
        std::vector<ldfReader::AcdDigi::AcdPmt>::const_iterator curReadout;
        if (readoutCol.size() > 2) log << MSG::DEBUG << "Too many readouts associated with this tile " << tileName << endreq;

        unsigned short pha[2] = {0,0};
        bool hitMapArr[2] = {false, false};
        //  Assumed passed low threshold - it looks like zero suppression is always on now
        bool acceptMapArr[2] = {true, true};
        bool cnoArr[2] = {false, false};
        Event::AcdDigi::Range range[2] = {Event::AcdDigi::LOW, Event::AcdDigi::LOW};
        Event::AcdDigi::ParityError error[4] = { Event::AcdDigi::NOERROR, Event::AcdDigi::NOERROR, Event::AcdDigi::NOERROR, Event::AcdDigi::NOERROR};

        for (curReadout = readoutCol.begin(); curReadout != readoutCol.end(); curReadout++) {
            int index = (curReadout->getSide() == ldfReader::AcdDigi::A) ? 0 : 1;
            pha[index] = (unsigned short) curReadout->getPha();
            hitMapArr[index] = curReadout->getHit();
            acceptMapArr[index] = curReadout->getAccept();

            range[index] = (curReadout->getRange() == 0) ? Event::AcdDigi::LOW : Event::AcdDigi::HIGH;
            error[index] = (curReadout->getParityError() == ldfReader::AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
            // set header Parity bit each PMT..since they're on different FREE boards
            error[index+2] = (curReadout->getHeaderParity() == ldfReader::AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
        }
        Event::AcdDigi *newDigi = new Event::AcdDigi(
            identsId, identsId.volId(), 0.0, pha, hitMapArr, acceptMapArr, cnoArr);


        newDigi->initLdfParameters(tileName, tileNumber, range, error);

        digiCol->push_back(newDigi);
        
    }

    return StatusCode::SUCCESS;
}

StatusCode LdfAcdDigiCnv::updateObj(int* , Event::AcdDigi* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfAcdDigiCnv");
    log << MSG::DEBUG << "LdfAcdDigiCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}

void LdfAcdDigiCnv::base10ToAcdId(unsigned int val, short &lay, short &face, 
                          short &row, short &col) {
    lay = val / 1000;
    val -= lay*1000;
    face = val / 100;
    val -= face*100;
    row = val / 10;
    val -= row*10;
    col = val;
}
