// File and Version Information:
//      $Header$
//
// Description:
//      LdfTkrDigiCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfTkrDigiCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "ldfReader/data/LatData.h"
#include "Event/Digi/TkrDigi.h"
#include "LdfBaseCnv.h"


class  LdfTkrDigiCnv : public LdfBaseCnv
//public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfTkrDigiCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfTkrDigiCnv(ISvcLocator* svc);

    virtual ~LdfTkrDigiCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_TkrDigi;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_TkrDigi;}

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
DECLARE_CONVERTER_FACTORY ( LdfTkrDigiCnv );

LdfTkrDigiCnv::LdfTkrDigiCnv(ISvcLocator* svc)  : LdfBaseCnv(classID(), svc)
//: Converter(TEST_StorageType, CLID_TkrDigi, svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/Digi/TkrDigiCol", objType(), "PASS");
    m_path = "/Event/Digi/TkrDigiCol";
}
/*
StatusCode LdfTkrDigiCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode LdfTkrDigiCnv::finalize() 
{
    return Converter::finalize();
}
*/

StatusCode LdfTkrDigiCnv::createObj(IOpaqueAddress* , 
                                    DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    //refpObject = new Event::EventHeader();
    Event::TkrDigiCol *digiCol = new Event::TkrDigiCol;
    refpObject = digiCol;
    
    // Retrieve the LAT data for this event and print out its contents
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    unsigned int iTower;
    for (iTower = 0; iTower < 16; iTower++) {
        ldfReader::TowerData *tem = myLatData->getTower(iTower);
        if (!tem) {
            //printf("Tem does not exist\n");
        } else {
            unsigned int i;
            for (i = 0; i < 18; i++) {
                ldfReader::TkrLayer *layer = tem->getTkrLayer(i);
                if (layer) {
                    
                    idents::TowerId tower = idents::TowerId(iTower);
                    int layerNum = i;
                    int ToT[2] = {0,0};
                    
                    Event::TkrDigi *xTdsLayer = 0, *yTdsLayer = 0;
                    
                    if (layer->getNumXstrips() > 0) {
                        ToT[0] = layer->getX_ToT(0);
                        ToT[1] = layer->getX_ToT(1);
                        xTdsLayer = new Event::TkrDigi(layerNum, idents::GlastAxis::X, tower, ToT);
                        std::vector<unsigned short> stripC0Col = layer->getXStripC0Col();
                        std::sort(stripC0Col.begin(), stripC0Col.end());
                        std::vector<unsigned short>::const_iterator stripC0;
                        for (stripC0 = stripC0Col.begin(); stripC0 != stripC0Col.end(); stripC0++) {
                            xTdsLayer->addC0Hit(*stripC0);
                        }
                        
                        std::vector<unsigned short> stripC1Col = layer->getXStripC1Col();
                        std::sort(stripC1Col.begin(), stripC1Col.end());
                        std::vector<unsigned short>::const_iterator stripC1;
                        for (stripC1 = stripC1Col.begin(); stripC1 != stripC1Col.end(); stripC1++) {
                            xTdsLayer->addC1Hit(*stripC1);
                        }

                        digiCol->push_back(xTdsLayer);
                    }
                    
                    if (layer->getNumYstrips() > 0) {
                        ToT[0] = layer->getY_ToT(0);
                        ToT[1] = layer->getY_ToT(1);
                        yTdsLayer = new Event::TkrDigi(layerNum, idents::GlastAxis::Y, tower, ToT);
                        std::vector<unsigned short> stripC0Col = layer->getYStripC0Col();
                        std::sort(stripC0Col.begin(), stripC0Col.end());
                        std::vector<unsigned short>::const_iterator stripC0;
                        for (stripC0 = stripC0Col.begin(); stripC0 != stripC0Col.end(); stripC0++) {
                            yTdsLayer->addC0Hit(*stripC0);
                        }
                        
                        std::vector<unsigned short> stripC1Col = layer->getYStripC1Col();
                        std::sort(stripC1Col.begin(), stripC1Col.end());
                        std::vector<unsigned short>::const_iterator stripC1;
                        for (stripC1 = stripC1Col.begin(); stripC1 != stripC1Col.end(); stripC1++) {
                            yTdsLayer->addC1Hit(*stripC1);
                        }

                        digiCol->push_back(yTdsLayer);
                    }                  
                }
            }           
        }
    }

    // sort by tower, layer, view
    std::sort(digiCol->begin(), digiCol->end(), Event::TkrDigi::digiLess());


    return StatusCode::SUCCESS;
}

//StatusCode LdfTkrDigiCnv::updateObj(int* , Event::TkrDigi* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

  //  MsgStream log(msgSvc(), "LdfTkrDigiCnv");
 //   log << MSG::DEBUG << "LdfTkrDigiCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


