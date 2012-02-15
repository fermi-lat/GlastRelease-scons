/* implementation file for GlastRandomSvc

gets adresses 
 of CLHEP random number engines used in Gleam shared libraries
 and sets seeds for them based on run and particle sequence
 number obtained from the MCHeader

 $Header$

 Author: Toby Burnett, Karl Young
*/

#include "GlastRandomObs.h"
#include "facilities/Util.h"

#include <iterator>
#include <fstream>
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataStoreItem.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/SmartDataPtr.h"


#include "CLHEP/Random/RandomEngine.h"

// for list of random engines
#include "CLHEP/Random/DualRand.h"
#include "CLHEP/Random/TripleRand.h"
#include "CLHEP/Random/RandEngine.h"
#include "CLHEP/Random/DRand48Engine.h"
#include "CLHEP/Random/Hurd160Engine.h"
#include "CLHEP/Random/Hurd288Engine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/RanshiEngine.h"
#include "CLHEP/Random/RandGauss.h"


/** @class RanGenFactoryBase
 *
* @brief abstract base class for list of HepRandomEngine factories 
*/
class RanGenFactoryBase {
public:
  virtual ~RanGenFactoryBase(){}
    virtual CLHEP::HepRandomEngine * create() =0;
    virtual std::string name()const =0;
};

/** @class RanGenFactory
 *
* @brief Templated class  to create a particular HepRandomEngine factory 
*/
template <class T> class RanGenFactory : public RanGenFactoryBase{ 
public:
    RanGenFactory(std::string name):m_name(name){}
    virtual CLHEP::HepRandomEngine * create(){return new T;}
    std::string name()const { return m_name;}
private:
    std::string m_name;

};

GlastRandomObs::GlastRandomObs():IToolSvc::Observer(),m_randomEngine("TripleRand"),
                               m_runNumber(-1),m_sequenceNumber(0),
                               m_autoSeed(false), m_engineNum(0)
{
    m_setFlagPointers.clear();
}

CLHEP::HepRandomEngine* GlastRandomObs::createEngine(std::string  engineName) 
{
    // Purpose and Method: Given a string containing a valid CLHEP
    //                     random engine type returns an instance 
    //
    // Inputs:  Name of CLHEP random engine to set
    // Outputs: pointer to a HepRandomEngine, or zero
    // Dependencies: None
    // Restrictions and Caveats:  None

    static RanGenFactoryBase*  factories[] = { 
        new RanGenFactory<CLHEP::DRand48Engine>("DRand48Engine")   ,
        new RanGenFactory<CLHEP::DualRand>("DualRand")             ,
        new RanGenFactory<CLHEP::TripleRand>("TripleRand")         ,
        new RanGenFactory<CLHEP::RandEngine>("RandEngine")         ,
        new RanGenFactory<CLHEP::HepJamesRandom>("HepJamesRandom") ,
        new RanGenFactory<CLHEP::HepJamesRandom>("JamesRandom")    ,
        new RanGenFactory<CLHEP::Hurd160Engine>("Hurd160Engine")   ,
        new RanGenFactory<CLHEP::Hurd288Engine>("Hurd288Engine")   ,
        new RanGenFactory<CLHEP::MTwistEngine>("MTwistEngine")     ,
        new RanGenFactory<CLHEP::RanecuEngine>("RanecuEngine")     ,
        new RanGenFactory<CLHEP::RanshiEngine>("RanshiEngine")     ,
        new RanGenFactory<CLHEP::RanluxEngine>("RanluxEngine")     ,
        new RanGenFactory<CLHEP::Ranlux64Engine>("Ranlux64Engine") 
    };

    for(unsigned int i=0; i< sizeof(factories)/sizeof(void*); ++i){
        if( engineName==factories[i]->name()) return factories[i]->create();
    }
    return 0;
}

//void GlastRandomObs::onRetrieve(const IAlgTool* tool) {
//	StatusCode status;//
//}

void GlastRandomObs::onCreate(const IAlgTool* tool) {

    //MsgStream log (msgSvc(), name() ); 
    IRandomAccess* ranacc ;
    StatusCode status =const_cast<IAlgTool*>(tool)->queryInterface(IRandomAccess::interfaceID(), (void**)&ranacc);
    if( status.isSuccess() )
    {   
        std::string tooltype = tool->type();
        // Set the Random engine by name
        CLHEP::HepRandomEngine* hr = createEngine(m_randomEngine);
        if( hr==0) 
        {
  //          log << MSG::ERROR << "Could not set up the random engine " 
            std::cout << "GlastRandomObs Could not set up the random engine "
                << m_randomEngine << std::endl;
//            return StatusCode::FAILURE;
        }
        //log << MSG::INFO << "Setting CLHEP Engine "<< m_randomEngine
        std::cout << "GlastRandomObs INFO Setting CLHEP Engine " 
                  << m_randomEngine
                  << " for " << tooltype << " at " << hr << std::endl;
        ++m_engineNum;
        CLHEP::HepRandomEngine* old = ranacc->setTheEngine(hr);
        for (EngineMap::iterator eit = m_engineMap.begin(); eit != m_engineMap.end(); ++eit) {
            if ( eit->second != old) continue;
            // MWK: I downgrade this from WARNING to INFO.
            // We know that linux/gcc instantiates only one engine,
            // and the fact is handled in the code.
           //log << MSG::INFO << " Previous engine ("<< old << ") for " 
           std::cout << "GlastRandomObs INFO Previous engine (" << old <<
              ") for " << tooltype << " was also set for "<< eit->first << std::endl;
           --m_engineNum;
           break;
         }

        m_engineMap[tooltype] = hr;  // store the RandomNum generator in the map

        IRandomAccess::SetFlag randset = ranacc->getRandSet();
        //log << MSG::INFO << "RandGauss setFlag at " << randset << endreq;
        std::cout << "RandGauss setFlag at " << randset << std::endl;
        m_setFlagPointers.push_back(randset);

   }
}



GlastRandomObs::~GlastRandomObs () 
{
    // Purpose and Method: finish up, delete the list of engines
    //
    // Inputs:  none
    // Outputs: none
    // Dependencies: None
    // Restrictions and Caveats:  None

    for (EngineMap::iterator dllEngine = m_engineMap.begin(); 
        dllEngine != m_engineMap.end(); ++dllEngine ) {
        delete dllEngine->second;       
    }
    m_engineMap.clear();
    m_setFlagPointers.clear();
}


