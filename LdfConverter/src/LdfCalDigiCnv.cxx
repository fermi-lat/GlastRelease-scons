// File and Version Information:
//      $Header$
//
// Description:
//      LdfCalDigiCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfCalDigiCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "Event/Digi/CalDigi.h"

#include "ldfReader/data/LatData.h"


class  LdfCalDigiCnv : public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfCalDigiCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfCalDigiCnv(ISvcLocator* svc);

    virtual ~LdfCalDigiCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_CalDigi;}
    static const unsigned char storageType() {return TEST_StorageType;}

    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_CalDigi;}

    /// Retrieve the class type of the data store the converter uses.
    // MSF: Masked to generate compiler error due to interface change
    virtual long repSvcType() const {return Converter::i_repSvcType();}

    /// Create the transient representation of an object.
    virtual StatusCode createObj(IOpaqueAddress* pAddress,DataObject*& refpObject);

    /// Methods to set and return the path in TDS for output of this converter
  //  virtual void setPath(const std::string& path) {m_path = path;}
    virtual const std::string& getPath() const    {return m_path;}

private:

    std::string m_path;

};

// Instantiation of a static factory class used by clients to create
// instances of this service
//static CnvFactory<MCEventCnv> s_factory;
//const ICnvFactory& MCEventCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY ( LdfCalDigiCnv );

LdfCalDigiCnv::LdfCalDigiCnv(ISvcLocator* svc)
: Converter(TEST_StorageType, CLID_CalDigi, svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    m_path = "/Event/Digi/CalDigiCol";
}

StatusCode LdfCalDigiCnv::initialize() {
    return Converter::initialize();
}

StatusCode LdfCalDigiCnv::finalize() {
    return Converter::finalize();
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
                    if (!readout) {
                        printf("LdfCalDigiCnv:  Faen! In BESTRANGE mode but first readout does not exist\n");
                        exit(1);
                    }
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
                        if (!readout) {
                            //printf("LdfCalDigiCnv:  Faen!  In ALLRANGE mode but not all range readouts are populated!\n");
                            //exit(1);
                            //I&T requests that we continue to process this data
                            // even if ranges are missing
                            continue;
                        }
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

//StatusCode LdfCalDigiCnv::updateObj(int* , Event::CalDigi* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.
//
//    MsgStream log(msgSvc(), "LdfCalDigiCnv");
//    log << MSG::DEBUG << "LdfCalDigiCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


