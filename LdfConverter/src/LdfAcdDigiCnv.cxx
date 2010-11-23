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
#include "LdfEvent/Gem.h"

#include "idents/AcdId.h"
#include "facilities/Util.h"

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

bool LdfAcdDigiCnv::createAcdDigi(Event::AcdDigiCol*& digiCol, unsigned int bitWord, int face, int offset, int maxRow, int maxCol) {
    bool created = false;
    int layer = 0;
    unsigned int iRow, iCol;
    for (iRow = 0; iRow <= maxRow; iRow++) {
        for (iCol = 0; iCol <= maxCol; iCol++) {
            // Right shift as necessary and check the bit
            if ( !((bitWord >> (iCol+offset*iRow)) & 1) ) continue;
            created = true;
            idents::AcdId identsId(layer, face, iRow, iCol);
            idents::VolumeIdentifier volId = identsId.volId();
            Event::AcdDigi *digi = new Event::AcdDigi(identsId, volId);
            digi->setNinja(true);
            digi->setGem(true);
            digiCol->push_back(digi);
        }
    }
    return created;
}

StatusCode LdfAcdDigiCnv::createObj(IOpaqueAddress* , DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    MsgStream log(msgSvc(), "LdfAcdDigiCnv");
    log << MSG::DEBUG << "Inside createObj" << endreq;
    std::map<idents::AcdId, Event::AcdDigi*> digiMapCol;

    Event::AcdDigiCol *digiCol = new Event::AcdDigiCol;
    refpObject = digiCol;

    // Retrieve the LAT data for this event and grab the ACD digi data
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    if (!myLatData) return StatusCode::SUCCESS;
 
    ldfReader::GemData ldfGem = myLatData->getGem();
    unsigned short cnoVector = ldfGem.cnoVector();
    /// This cnoVector is described in 
    /// http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/GEM.pdf 
    /// in Figure 124 and read the description in Section 1.6.5

    // AKA GEM hitmap.  See the GEM documetation cited above Section 4.9
    const ldfReader::GemDataTileList& tileList = ldfGem.tileList();
    unsigned short xzm = tileList.XZM(); // -Y side 200s
    unsigned short xzp = tileList.XZP(); // +Y side 400s
    unsigned short yzm = tileList.YZM(); // -X side 100s
    unsigned short yzp = tileList.YZP(); // +X side 300s
    unsigned xy = tileList.XY();         // top
    unsigned short rbn = tileList.RBN(); // ribons
    unsigned short na = tileList.NA();   // Not connected

    // Using the GEM map, create AcdDigi objects
    int face=0, offset = 5;
    createAcdDigi(digiCol, xy, face, offset, 4, 4);

    // sides
    face = 1;
    createAcdDigi(digiCol, yzm, face, offset, 3, 4);

    face = 2;
    createAcdDigi(digiCol, xzm, face, offset, 3, 4);

    face = 3;
    createAcdDigi(digiCol, yzp, face, offset, 3, 4);

    face = 4;
    createAcdDigi(digiCol, xzp, face, offset, 3, 4); 

    face = 5;
    offset = 4;
    createAcdDigi(digiCol, rbn, face, offset, 0, 3);

    face = 6;
    unsigned short rbn2 = rbn >> 4;
    createAcdDigi(digiCol, rbn2, face, offset, 0, 3);

    // handle NAs separately
    int layer = 1;
    face = 0;
    int row = 0, col = 0;
    int iBit;
    for (iBit = 0; iBit <= 10; iBit++) {
        if (!(na & (1 << iBit))) continue;
        idents::VolumeIdentifier volId;
        if (iBit == 10) {
            row = 1;
            col = 0;
        } else {
            col = iBit;
        }

        idents::AcdId acdId(layer, face, row, col);
        
        std::string tileName = "NA";
        std::string numStr;
        facilities::Util::itoa(iBit,numStr);
        tileName += numStr;
    }

    // Create a map for easy retrieval of the created AcdDigis
    std::vector<Event::AcdDigi*>::const_iterator vecIt;
    for (vecIt = digiCol->begin(); vecIt != digiCol->end(); vecIt++){
        digiMapCol[(*vecIt)->getId()] = *vecIt;
    }

    const std::map<const char*, ldfReader::AcdDigi*> acdCol = myLatData->getAcdCol();
    if (digiCol->size() != acdCol.size()) {
        log << MSG::DEBUG << "GEM created AcdDigis: " << digiCol->size()
            << " AcdDigis found in hitmap: " << acdCol.size() << endreq;
    }
    std::map<const char*, ldfReader::AcdDigi*>::const_iterator thisAcdDigi;
    
    log << MSG::DEBUG << "found acddigis: " << acdCol.size() << endreq;

    // Now loop over the original AcdDigis and search for the AcdDigi in the digiCol
    for (thisAcdDigi = acdCol.begin(); thisAcdDigi != acdCol.end(); thisAcdDigi++) {
        const char *tileName = thisAcdDigi->second->getTileName();
        int tileNumber = thisAcdDigi->second->getTileNumber();
        unsigned int tileId = thisAcdDigi->second->getTileId();
        log << MSG::DEBUG << "tileId " << tileId << endreq;
        short layer, face, row, col;
        base10ToAcdId(tileId, layer, face, row, col);
        idents::AcdId identsId(layer, face, row, col);
        // avoid volId creation for NA
        idents::VolumeIdentifier volId;
        // Test for N/A
        if ((tileId == 0) && (strcmp(tileName, "000") != 0) ) 
            identsId.na(1);
        if (!identsId.na())
            volId = identsId.volId();

        log << MSG::DEBUG << "idents id " << identsId.id() << endreq; 

        const std::vector<ldfReader::AcdDigi::AcdPmt> readoutCol = thisAcdDigi->second->getReadout();
        log << MSG::DEBUG << "got readout col size " << readoutCol.size() << endreq;
        std::vector<ldfReader::AcdDigi::AcdPmt>::const_iterator curReadout;
        if (readoutCol.size() > 2) log << MSG::DEBUG << "Too many readouts associated with this tile " << tileName << endreq;

        unsigned short pha[2] = {0,0};
        bool hitMapArr[2] = {false, false};
        bool acceptMapArr[2] = {false, false};
        bool cnoArr[2] = {false, false};
        Event::AcdDigi::Range range[2] = {Event::AcdDigi::LOW, Event::AcdDigi::LOW};
        Event::AcdDigi::ParityError error[4] = { Event::AcdDigi::NOERROR, Event::AcdDigi::NOERROR, Event::AcdDigi::NOERROR, Event::AcdDigi::NOERROR};

        for (curReadout = readoutCol.begin(); curReadout != readoutCol.end(); curReadout++) {
            int index = (curReadout->getSide() == ldfReader::AcdDigi::A) ? 0 : 1;
            pha[index] = (unsigned short) curReadout->getPha();
            hitMapArr[index] = curReadout->getHit();
            acceptMapArr[index] = curReadout->getAccept();
            cnoArr[index] = ((cnoVector >> curReadout->getCable()) & 1) ? true : false;

            range[index] = (curReadout->getRange() == 0) ? Event::AcdDigi::LOW : Event::AcdDigi::HIGH;
            error[index] = (curReadout->getParityError() == ldfReader::AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
            // set header Parity bit each PMT..since they're on different FREE boards
            error[index+2] = (curReadout->getHeaderParity() == ldfReader::AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
        }

        std::map<idents::AcdId, Event::AcdDigi*>::iterator findDigi = digiMapCol.find(identsId);
        if (findDigi == digiMapCol.end()) {

            log << MSG::DEBUG << "Creating digi obj " << endreq;
            Event::AcdDigi *newDigi = new Event::AcdDigi (
                identsId, volId, 0.0, pha, hitMapArr, acceptMapArr, cnoArr);


            newDigi->initLdfParameters(tileName, tileNumber, range, error);

            digiCol->push_back(newDigi);
        } else {
            findDigi->second->init(pha, hitMapArr, acceptMapArr, cnoArr);
            findDigi->second->initLdfParameters(tileName, tileNumber, range, error);
            if ((hitMapArr[0] == true) && (acceptMapArr[0] == false) )
                findDigi->second->setNinja(true);
            else if ((hitMapArr[1] == true) && (acceptMapArr[1] == false) )
                findDigi->second->setNinja(true);
            else
                findDigi->second->setNinja(false);
        }
        
    }

    log << MSG::DEBUG << "Returning" << endreq;
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
