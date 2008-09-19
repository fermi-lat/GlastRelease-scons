#include <LdfEvent/LsfMetaEvent.h>
#include "LdfEvent/LsfCcsds.h"
#include "OnboardFilterTds/ObfFilterStatus.h"

#include <digiRootData/GemTime.h>
#include <digiRootData/TimeTone.h>
#include <digiRootData/LsfTime.h>
#include <digiRootData/RunInfo.h>
#include <digiRootData/DatagramInfo.h>
#include <digiRootData/GemScalers.h>
#include <digiRootData/Configuration.h>
#include <digiRootData/MetaEvent.h>
#include <digiRootData/Ccsds.h>

#include <RootConvert/Utilities/Toolkit.h>

#include "TObjArray.h"
#include "TCollection.h"


namespace RootPersistence {

  void convert( const lsfData::GemTime& tdsObj, GemTime& rootObj) {
    rootObj.initialize(tdsObj.hacks(),tdsObj.ticks());
  }; 
  void convert( const GemTime& rootObj, lsfData::GemTime& tdsObj) {    
    tdsObj.set(rootObj.hacks(),rootObj.ticks());
  };
  
  void convert( const lsfData::TimeTone& tdsObj, TimeTone& rootObj) {
    GemTime rootHack;
    convert(tdsObj.timeHack(),rootHack);
    rootObj.initialize(tdsObj.incomplete(),tdsObj.timeSecs(),tdsObj.flywheeling(),tdsObj.flags(),rootHack);
  };
  void convert( const TimeTone& rootObj, lsfData::TimeTone& tdsObj) {
    lsfData::GemTime tdsHack;
    convert(rootObj.timeHack(),tdsHack);
    tdsObj.set(rootObj.incomplete(),rootObj.timeSecs(),rootObj.flywheeling(),rootObj.flags(),tdsHack);
  };
  
  void convert( const lsfData::Time& tdsObj, LsfTime& rootObj) {
    TimeTone current;
    TimeTone previous;
    GemTime hack;
    convert(tdsObj.current(),current);
    convert(tdsObj.previous(),previous);
    convert(tdsObj.timeHack(),hack);
    rootObj.initialize(current,previous,hack,tdsObj.timeTicks());    
  };
  void convert( const LsfTime& rootObj, lsfData::Time& tdsObj) {
    lsfData::TimeTone current;
    lsfData::TimeTone previous;
    lsfData::GemTime hack;
    convert(rootObj.current(),current);
    convert(rootObj.previous(),previous);
    convert(rootObj.timeHack(),hack);
    tdsObj.set(current,previous,hack,rootObj.timeTicks());
  };

  void convert( const lsfData::RunInfo& tdsObj, RunInfo& rootObj) {
    rootObj.initialize(tdsObj.platform(),tdsObj.dataOrigin(),tdsObj.id(),tdsObj.startTime(),tdsObj.dataTransferId());
  };

  void convert( const RunInfo& rootObj, lsfData::RunInfo& tdsObj) {
    tdsObj.set(rootObj.platform(),rootObj.dataOrigin(),rootObj.id(),rootObj.startTime(),rootObj.dataTransferId());
  };
   
  void convert( const lsfData::DatagramInfo& tdsObj, DatagramInfo& rootObj) {
    rootObj.initialize(tdsObj.openAction(),tdsObj.openReason(),
		       tdsObj.crate(),tdsObj.mode(),
		       tdsObj.closeAction(),tdsObj.closeReason(),
		       tdsObj.datagrams(),tdsObj.modeChanges());    
  };
  void convert( const DatagramInfo& rootObj, lsfData::DatagramInfo& tdsObj) {
    tdsObj.set(rootObj.openAction(),rootObj.openReason(),
	       rootObj.crate(),rootObj.mode(),
	       rootObj.closeAction(),rootObj.closeReason(),
	       rootObj.datagrams(),rootObj.modeChanges());
  };

  void convert( const lsfData::GemScalers& tdsObj, GemScalers& rootObj) {
    rootObj.initialize(tdsObj.elapsed(),tdsObj.livetime(),
		       tdsObj.prescaled(),tdsObj.discarded(),
		       tdsObj.sequence(),tdsObj.deadzone());
  };
  void convert( const GemScalers& rootObj, lsfData::GemScalers& tdsObj) {
    tdsObj.set(rootObj.elapsed(),rootObj.livetime(),
	       rootObj.prescaled(),rootObj.discarded(),
	       rootObj.sequence(),rootObj.deadzone());
  };

  // this handles the polymorphism of Configuration
  void convert( const lsfData::Configuration* tdsPtr, Configuration*& rootPtr ) {
    rootPtr = 0;
    if ( tdsPtr == 0 ) return;
    const lsfData::LpaConfiguration* asLpaTds(0);
    switch ( tdsPtr->type() ) {
    case enums::Lsf::LPA:
      asLpaTds = tdsPtr->castToLpaConfig();
      rootPtr = new LpaConfiguration(asLpaTds->hardwareKey(),asLpaTds->softwareKey());
      break;
    default:
      return;
    }
  };

  void convert( const Configuration* rootPtr, lsfData::Configuration*& tdsPtr) {
    tdsPtr = 0;
    if ( rootPtr == 0 ) return;
    switch ( rootPtr->runType() ) {
    case enums::Lsf::LPA:
      break;
    default:
      return;
    }
  };

  void convert( const lsfData::LpaConfiguration* tdsObj,
                LpaConfiguration& rootObj) {

      rootObj.initialize(tdsObj->hardwareKey(), tdsObj->softwareKey());
  };

  void convert( const LpaConfiguration* rootObj,
                lsfData::LpaConfiguration*& tdsObj) {

      tdsObj->set(rootObj->hardwareKey(), rootObj->softwareKey());
  };

  void convert( const lsfData::LciAcdConfiguration* tdsObj, 
                LciAcdConfiguration& rootObj) {
      lsfData::LciAcdConfiguration::AcdTrigger triggerTds = tdsObj->trigger();

      LciAcdConfiguration::AcdTrigger triggerRoot(triggerTds.veto(),
         triggerTds.vetoVernier(), triggerTds.highDiscrim());

      lsfData::Channel channelTds = tdsObj->channel();

      Channel channelRoot(channelTds.single(), channelTds.all(), 
          channelTds.latc());

      rootObj.initialize(tdsObj->injected(), tdsObj->threshold(),
          tdsObj->biasDac(), tdsObj->holdDelay(), tdsObj->hitmapDelay(),
          tdsObj->range(), triggerRoot, channelRoot);

      rootObj.setSoftwareKey(tdsObj->softwareKey());
      rootObj.setWriteCfg(tdsObj->writeCfg());
      rootObj.setReadCfg(tdsObj->readCfg());
      rootObj.setPeriod(tdsObj->period());
      rootObj.setFlags(tdsObj->flags());

  };

  void convert( const LciAcdConfiguration* rootObj,
                lsfData::LciAcdConfiguration*& tdsObj) {

      LciAcdConfiguration::AcdTrigger triggerRoot = rootObj->trigger();
      lsfData::LciAcdConfiguration::AcdTrigger triggerTds(triggerRoot.veto(),
          triggerRoot.vetoVernier(), triggerRoot.highDiscrim());

      Channel chRoot = rootObj->channel();
      lsfData::Channel chTds(chRoot.single(), chRoot.all(), chRoot.latc());

      tdsObj->set(rootObj->injected(), rootObj->threshold(),
         rootObj->biasDac(), rootObj->holdDelay(), rootObj->hitmapDelay(),
         rootObj->range(), triggerTds, chTds);
      tdsObj->setSoftwareKey(rootObj->softwareKey());
      tdsObj->setWriteCfg(rootObj->writeCfg());
      tdsObj->setReadCfg(rootObj->readCfg());
      tdsObj->setPeriod(rootObj->period());
      tdsObj->setFlags(rootObj->flags());
  };

  void convert( const lsfData::LciCalConfiguration* tdsObj,
                LciCalConfiguration& rootObj) {
      lsfData::LciCalConfiguration::CalTrigger trigTds = tdsObj->trigger();

      LciCalConfiguration::CalTrigger trigRoot(trigTds.le(), 
                     trigTds.lowTrgEna(), trigTds.he(), trigTds.highTrgEna());

      lsfData::Channel chTds = tdsObj->channel();
      Channel chRoot(chTds.single(), chTds.all(), chTds.latc());

      rootObj.initialize(tdsObj->uld(), tdsObj->injected(), tdsObj->delay(),
              tdsObj->firstRange(), tdsObj->threshold(),tdsObj->calibGain(),
              tdsObj->highCalEna(),tdsObj->highRngEna(),tdsObj->highGain(),
              tdsObj->lowCalEna(),tdsObj->lowRngEna(),tdsObj->lowGain(),
              trigRoot,chRoot);
      rootObj.setSoftwareKey(tdsObj->softwareKey());
      rootObj.setWriteCfg(tdsObj->writeCfg());
      rootObj.setReadCfg(tdsObj->readCfg());
      rootObj.setPeriod(tdsObj->period());
      rootObj.setFlags(tdsObj->flags());
  };

  void convert( const LciCalConfiguration* rootObj,
                lsfData::LciCalConfiguration*& tdsObj) {

      LciCalConfiguration::CalTrigger trigRoot = rootObj->trigger();
      lsfData::LciCalConfiguration::CalTrigger trigTds(trigRoot.le(),
          trigRoot.lowTrgEna(), trigRoot.he(), trigRoot.highTrgEna());
      Channel chRoot = rootObj->channel();
      lsfData::Channel chTds(chRoot.single(), chRoot.all(), chRoot.latc());
      tdsObj->set(rootObj->uld(), rootObj->injected(), rootObj->delay(),
              rootObj->firstRange(),
              rootObj->threshold(),rootObj->calibGain(),rootObj->highCalEna(),
              rootObj->highRngEna(),rootObj->highGain(),rootObj->lowCalEna(),
              rootObj->lowRngEna(),rootObj->lowGain(),trigTds,chTds);
      tdsObj->setSoftwareKey(rootObj->softwareKey());
      tdsObj->setWriteCfg(rootObj->writeCfg());
      tdsObj->setReadCfg(rootObj->readCfg());
      tdsObj->setPeriod(rootObj->period());
      tdsObj->setFlags(rootObj->flags());
  };

  void convert( const lsfData::LciTkrConfiguration* tdsObj,
                LciTkrConfiguration& rootObj) {

      lsfData::Channel chTds = tdsObj->channel();
      Channel chRoot(chTds.single(), chTds.all(), chTds.latc());
      rootObj.initialize(tdsObj->injected(),tdsObj->delay(),
              tdsObj->threshold(), tdsObj->splitLow(), tdsObj->splitHigh(),
              chRoot);
      rootObj.setSoftwareKey(tdsObj->softwareKey());
      rootObj.setWriteCfg(tdsObj->writeCfg());
      rootObj.setReadCfg(tdsObj->readCfg());
      rootObj.setPeriod(tdsObj->period());
      rootObj.setFlags(tdsObj->flags());
  };

  void convert( const LciTkrConfiguration* rootObj,
                lsfData::LciTkrConfiguration*& tdsObj) {

      Channel chRoot = rootObj->channel();
      lsfData::Channel chTds(chRoot.single(), chRoot.all(), chRoot.latc());

      tdsObj->set(rootObj->injected(),rootObj->delay(),
              rootObj->threshold(), rootObj->splitLow(),rootObj->splitHigh(),
              chTds);
      tdsObj->setSoftwareKey(rootObj->softwareKey());
      tdsObj->setWriteCfg(rootObj->writeCfg());
      tdsObj->setReadCfg(rootObj->readCfg());
      tdsObj->setPeriod(rootObj->period());
      tdsObj->setFlags(rootObj->flags());
  };

  void convert( const LpaKeys* rootObj,
                lsfData::LpaKeys*& tdsObj){

    
    tdsObj->setSbs(rootObj->sbs());
    tdsObj->setLpa_db(rootObj->lpa_db());
    tdsObj->setLATC_master(rootObj->LATC_master());
    tdsObj->setLATC_ignore(rootObj->LATC_ignore());

  };

  void convert( const lsfData::LpaKeys* tdsObj,
                LpaKeys& rootObj) {

    rootObj.setSbs(tdsObj->sbs());
    rootObj.setLpa_db(tdsObj->lpa_db());
    rootObj.setLATC_master(tdsObj->LATC_master());
    rootObj.setLATC_ignore(tdsObj->LATC_ignore());
  };

  void convert( const LciKeys* rootObj,
                lsfData::LciKeys*& tdsObj){
    tdsObj->setLCI_script(rootObj->LCI_script());
    tdsObj->setLATC_master(rootObj->LATC_master());
    tdsObj->setLATC_ignore(rootObj->LATC_ignore());
  };

  void convert( const lsfData::LciKeys* tdsObj,
                LciKeys& rootObj) {

    rootObj.setLCI_script(tdsObj->LCI_script());
    rootObj.setLATC_master(tdsObj->LATC_master());
    rootObj.setLATC_ignore(tdsObj->LATC_ignore());
  };

  void convert( const LsfEvent::MetaEvent& tdsObj, MetaEvent& rootObj) {
    RunInfo run;
    DatagramInfo datagram;
    GemScalers scalers;
    LsfTime time;
    Configuration* config(0);
    convert(tdsObj.run(),run);
    convert(tdsObj.datagram(),datagram);
    convert(tdsObj.scalers(),scalers);
    convert(tdsObj.time(),time);
    convert(tdsObj.configuration(),config);

    rootObj.initialize(run,datagram,scalers,time,*config);
    rootObj.setRunType(tdsObj.runType());

    if ( tdsObj.runType() == enums::Lsf::LPA ) {
      const lsfData::LpaConfiguration *lpa = tdsObj.configuration()->castToLpaConfig();
      LpaConfiguration lpaRoot(lpa->hardwareKey(), lpa->softwareKey());
      rootObj.setLpaConfiguration(lpaRoot);
    } else if (tdsObj.runType() == enums::Lsf::AcdLCI) {
      const lsfData::LciAcdConfiguration* acdTds = tdsObj.configuration()->castToLciAcdConfig();
      LciAcdConfiguration acdRoot;
      convert(acdTds, acdRoot);
      rootObj.setLciAcdConfiguration(acdRoot);
    } else if (tdsObj.runType() == enums::Lsf::CalLCI) {
      const lsfData::LciCalConfiguration *calTds = tdsObj.configuration()->castToLciCalConfig();
      LciCalConfiguration calRoot;
      convert(calTds, calRoot);
      rootObj.setLciCalConfiguration(calRoot);
    } else if (tdsObj.runType() == enums::Lsf::TkrLCI ) {
      const lsfData::LciTkrConfiguration *tkrTds = tdsObj.configuration()->castToLciTkrConfig();
      LciTkrConfiguration tkrRoot;
      convert(tkrTds, tkrRoot);
      rootObj.setLciTkrConfiguration(tkrRoot);
    }

    if (config) delete config;

    if (tdsObj.keyType() == enums::Lsf::LpaKeys) {
      const lsfData::LpaKeys *lpaKeysTds = tdsObj.keys()->castToLpaKeys();
      LpaKeys lpaKeysRoot;
      convert(lpaKeysTds, lpaKeysRoot);
      rootObj.setLpaKeys(lpaKeysRoot);
    } else if (tdsObj.keyType() == enums::Lsf::LciKeys) {
      const lsfData::LciKeys *lciKeysTds = tdsObj.keys()->castToLciKeys();
      LciKeys lciKeysRoot;
      convert(lciKeysTds, lciKeysRoot);
      rootObj.setLciKeys(lciKeysRoot);
    }

    rootObj.setMootKey(tdsObj.mootKey());
    rootObj.setMootAlias(tdsObj.mootAlias());

    rootObj.setCompressionLevel(tdsObj.compressionLevel());
    rootObj.setCompressedSize(tdsObj.compressedSize());

    if (tdsObj.passthruFilter()) {
        const lsfData::PassthruHandler* tdsPass = tdsObj.passthruFilter();
        LpaPassthruFilter *pass = new LpaPassthruFilter;
        pass->initialize(tdsPass->masterKey(),tdsPass->cfgKey(),
            tdsPass->cfgId(),tdsPass->state(),tdsPass->prescaler(),
            tdsPass->version(),tdsPass->id(),tdsPass->has(),
            tdsPass->prescaleFactor());
        if(tdsPass->rsd())
            pass->setStatusWord(tdsPass->rsd()->status());
       rootObj.addPassthru(pass);
    }
    if (tdsObj.gammaFilter()) {
       const lsfData::GammaHandler *tdsGam = tdsObj.gammaFilter();
       LpaGammaFilter *gam =0;
       switch(tdsGam->version()) {
         case 0:
             gam = new LpaGammaFilterV0;
             break;
         case 1:
             gam = new LpaGammaFilterV1;
             break;
         case 2:
             gam = new LpaGammaFilterV2;
             break;
         case 3:
             gam = new LpaGammaFilterV3;
             break;
         default:
             gam = 0;
       }
       //LpaGammaFilter *gam = new LpaGammaFilter;
       gam->initialize(tdsGam->masterKey(),tdsGam->cfgKey(),
            tdsGam->cfgId(),tdsGam->state(),tdsGam->prescaler(),
            tdsGam->version(),tdsGam->id(),tdsGam->has(),
            tdsGam->prescaleFactor());
       if (tdsGam->rsd())
           gam->set(tdsGam->rsd()->status(), tdsGam->rsd()->stage(), 
                     tdsGam->rsd()->energyValid(), tdsGam->rsd()->energyInLeus());
        rootObj.addGamma(gam);
    }
    if (tdsObj.mipFilter()){
       const lsfData::MipHandler *tdsMip = tdsObj.mipFilter();
       LpaMipFilter *mip = new LpaMipFilter;
       mip->initialize(tdsMip->masterKey(),tdsMip->cfgKey(),
                tdsMip->cfgId(),tdsMip->state(),tdsMip->prescaler(),
                tdsMip->version(),tdsMip->id(),tdsMip->has(),
                tdsMip->prescaleFactor());
       if (tdsMip->rsd())
           mip->setStatusWord(tdsMip->rsd()->status());
      rootObj.addMip(mip);
    }
    if (tdsObj.hipFilter()) {
       const lsfData::HipHandler *tdsHip = tdsObj.hipFilter();
       LpaHipFilter *hip = new LpaHipFilter;
       hip->initialize(tdsHip->masterKey(),tdsHip->cfgKey(),
             tdsHip->cfgId(),tdsHip->state(),tdsHip->prescaler(),
             tdsHip->version(),tdsHip->id(),tdsHip->has(),
             tdsHip->prescaleFactor());
       if (tdsHip->rsd())
           hip->setStatusWord(tdsHip->rsd()->status());
       rootObj.addHip(hip);
    }
    if (tdsObj.dgnFilter()) {
       const lsfData::DgnHandler *tdsDgn = tdsObj.dgnFilter();
       LpaDgnFilter *dgn = new LpaDgnFilter;
      
       dgn->initialize(tdsDgn->masterKey(),tdsDgn->cfgKey(),
             tdsDgn->cfgId(),tdsDgn->state(),tdsDgn->prescaler(),
             tdsDgn->version(),tdsDgn->id(),tdsDgn->has(),
             tdsDgn->prescaleFactor());
       if (tdsDgn->rsd())
           dgn->setStatusWord(tdsDgn->rsd()->status());
       rootObj.addDgn(dgn);
    }
/*  Holding off on this guy for now
    if(tdsObj.lpaHandler()) {
      const lsfData::LpaHandler *tdsHandler = tdsObj.lpaHandler();
      ILpaHandler* handler = new ILpaHandler;
      handler->initialize(tdsHandler->masterKey(),tdsHandler->cfgKey(),
               tdsHandler->cfgId(),tdsHandler->state(),
               tdsHandler->prescaler(), tdsHandler->version(),
               tdsHandler->id(),tdsHandler->has());
       rootObj.addLpaHandler(handler);
     }
*/


  };

  void convert( const MetaEvent& rootObj, LsfEvent::MetaEvent& tdsObj) {
    lsfData::RunInfo run;
    lsfData::DatagramInfo datagram;
    lsfData::GemScalers scalers;
    lsfData::Time time;
    lsfData::Configuration* config(0);
    lsfData::LsfKeys* keys(0);
    convert(rootObj.run(),run);
    convert(rootObj.datagram(),datagram);
    convert(rootObj.scalers(),scalers);
    convert(rootObj.time(),time);

    if ( rootObj.runType() == enums::Lsf::LPA ) {
        const LpaConfiguration* lpaRoot = rootObj.lpaConfiguration();
        lsfData::LpaConfiguration *lpaTds = new lsfData::LpaConfiguration;
        convert(lpaRoot,lpaTds);
        config = lpaTds->clone();
        delete lpaTds;
    } else if (rootObj.runType() == enums::Lsf::AcdLCI) {
        const LciAcdConfiguration* acdRoot = rootObj.lciAcdConfiguration();
        lsfData::LciAcdConfiguration *acdTds = new lsfData::LciAcdConfiguration;
        convert(acdRoot, acdTds);
        config = acdTds->clone();
        delete acdTds;
    } else if (rootObj.runType() == enums::Lsf::CalLCI) {
        const LciCalConfiguration* calRoot = rootObj.lciCalConfiguration();
        lsfData::LciCalConfiguration *calTds = new lsfData::LciCalConfiguration;
        convert(calRoot, calTds);
        config = calTds->clone();
        delete calTds;
    } else if (rootObj.runType() == enums::Lsf::TkrLCI ) {
        const LciTkrConfiguration* tkrRoot = rootObj.lciTkrConfiguration();
        lsfData::LciTkrConfiguration *tkrTds = new lsfData::LciTkrConfiguration;
        convert(tkrRoot, tkrTds);
        config = tkrTds->clone();
        delete tkrTds;
    }

    if( rootObj.keyType() == enums::Lsf::LpaKeys) {
      const LpaKeys *lpaKeysRoot = rootObj.lpaKeys();
      lsfData::LpaKeys *lpaKeysTds = new lsfData::LpaKeys;
      convert(lpaKeysRoot,lpaKeysTds);
      keys = lpaKeysTds;
    } else if (rootObj.keyType() == enums::Lsf::LciKeys) {
      const LciKeys *lciKeysRoot = rootObj.lciKeys();
      lsfData::LciKeys *lciKeysTds = new lsfData::LciKeys;
      convert(lciKeysRoot, lciKeysTds);
      keys = lciKeysTds;
    }

    tdsObj.set(run,datagram,scalers,time);
    if(config) {
        tdsObj.setConfiguration(*config);
        delete config;
    }
    if (keys) {
        tdsObj.setKeys(*keys);
        delete keys;
    }

    tdsObj.setMootKey(rootObj.mootKey());
    tdsObj.setMootAlias(rootObj.mootAlias());

    tdsObj.setCompressionLevel(rootObj.compressionLevel());
    tdsObj.setCompressedSize(rootObj.compressedSize());

    if (rootObj.lpaHandler().getPassthruFilter()) {
        const LpaPassthruFilter *handlerIt = rootObj.lpaHandler().getPassthruFilter();
        lsfData::PassthruHandler pass;
        pass.set(handlerIt->getMasterKey(),handlerIt->getCfgKey(),
                 handlerIt->getCfgId(), handlerIt->getState(),
                 handlerIt->getPrescaler(),handlerIt->getVersion(),
                 handlerIt->getId(),handlerIt->has(),
                 handlerIt->getPrescaleFactor());
        if (handlerIt->has())
            pass.setStatus(handlerIt->getStatusWord());
        tdsObj.addPassthruHandler(pass);
    }
    if (rootObj.lpaHandler().getGammaFilter()) {
        const LpaGammaFilter *handlerIt = rootObj.lpaHandler().getGammaFilter();
        lsfData::GammaHandler gam;
        gam.set(handlerIt->getMasterKey(),handlerIt->getCfgKey(),
                 handlerIt->getCfgId(), handlerIt->getState(),
                 handlerIt->getPrescaler(),handlerIt->getVersion(),
                 handlerIt->getId(),handlerIt->has(),
                 handlerIt->getPrescaleFactor());
        if (handlerIt->has())
            gam.setStatus(handlerIt->getStatusWord(),
                      handlerIt->getStage(), handlerIt->getEnergyValid(),
                      handlerIt->getEnergyInLeus());
        tdsObj.addGammaHandler(gam);  
    }
    if (rootObj.lpaHandler().getMipFilter()) {
        const LpaMipFilter *handlerIt = rootObj.lpaHandler().getMipFilter();
        lsfData::MipHandler mip;
        mip.set(handlerIt->getMasterKey(),handlerIt->getCfgKey(),
                 handlerIt->getCfgId(), handlerIt->getState(),
                 handlerIt->getPrescaler(),handlerIt->getVersion(),
                 handlerIt->getId(),handlerIt->has(),
                 handlerIt->getPrescaleFactor());
        if (handlerIt->has())
            mip.setStatus(handlerIt->getStatusWord());
        tdsObj.addMipHandler(mip);
    }
    if (rootObj.lpaHandler().getHipFilter()) {
        const LpaHipFilter *handlerIt = rootObj.lpaHandler().getHipFilter();
        lsfData::HipHandler hip;
        hip.set(handlerIt->getMasterKey(),handlerIt->getCfgKey(),
                 handlerIt->getCfgId(), handlerIt->getState(),
                 handlerIt->getPrescaler(),handlerIt->getVersion(),
                 handlerIt->getId(),handlerIt->has(),
                 handlerIt->getPrescaleFactor());
        if (handlerIt->has())
            hip.setStatus(handlerIt->getStatusWord());
        tdsObj.addHipHandler(hip);
     }
    if (rootObj.lpaHandler().getDgnFilter()) {
        const LpaDgnFilter *handlerIt = rootObj.lpaHandler().getDgnFilter();
        lsfData::DgnHandler dgn;
        dgn.set(handlerIt->getMasterKey(),handlerIt->getCfgKey(),
                 handlerIt->getCfgId(), handlerIt->getState(),
                 handlerIt->getPrescaler(),handlerIt->getVersion(),
                 handlerIt->getId(),handlerIt->has(),
                 handlerIt->getPrescaleFactor());
        if (handlerIt->has())
            dgn.setStatus(handlerIt->getStatusWord());
        tdsObj.addDgnHandler(dgn);
    }



  };

  void convert( const LsfEvent::LsfCcsds& tdsObj, Ccsds& rootObj) {
      rootObj.initialize(tdsObj.scid(), tdsObj.apid(), tdsObj.utc());
  };

  void convert (const Ccsds& rootObj, LsfEvent::LsfCcsds& tdsObj) {
      tdsObj.initialize(rootObj.getScid(), rootObj.getApid(), rootObj.getUtc());
  };


  void convert( const OnboardFilterTds::ObfFilterStatus& obf, MetaEvent& meta) {

    // Now step through and initialize the status objects one by one
    const OnboardFilterTds::IObfStatus* tdsStatus = 0;

    bool has = true;   // Gleam OBF always creates RSD
    unsigned int masterKey = LSF_INVALID_UINT; 
    unsigned int cfgKey = LSF_INVALID_UINT;
    unsigned int cfgId = LSF_INVALID_UINT;
    unsigned int version = 0;
    unsigned int prescaleFactor = LSF_INVALID_UINT;

    // We do this one by one explicitly for now. Start with the results of the gamma filter
    if ((tdsStatus = obf.getFilterStatus(OnboardFilterTds::ObfFilterStatus::GammaFilter)))
    {
        const OnboardFilterTds::ObfGammaStatus* gamObj = dynamic_cast<const OnboardFilterTds::ObfGammaStatus*>(tdsStatus);

       LpaGammaFilter *gam = new LpaGammaFilter;
       gam->initialize(masterKey,cfgKey,
            cfgId,enums::Lsf::RsdState(gamObj->getState()),
            enums::Lsf::LeakedPrescaler(gamObj->getPrescalerWord()),
            version,enums::Lsf::GAMMA,has, prescaleFactor);
       gam->set(gamObj->getStatusWord(), gamObj->getStage(), 
                     gamObj->getStage()&1, gamObj->getEnergyInLeus());
        meta.addGamma(gam);
    }

    // HFC Filter is next
    if ((tdsStatus = obf.getFilterStatus(OnboardFilterTds::ObfFilterStatus::HIPFilter)))
    {
        const OnboardFilterTds::ObfHipStatus* hipObj = dynamic_cast<const OnboardFilterTds::ObfHipStatus*>(tdsStatus);
       LpaHipFilter *hip = new LpaHipFilter;
       hip->initialize(masterKey,cfgKey,
            cfgId,enums::Lsf::RsdState(hipObj->getState()),
            enums::Lsf::LeakedPrescaler(hipObj->getPrescalerWord()),
            version,enums::Lsf::HIP,has, prescaleFactor);
       hip->setStatusWord(hipObj->getStatusWord());
       meta.addHip(hip);
    }

    // MIP Filter is next
    if ((tdsStatus = obf.getFilterStatus(OnboardFilterTds::ObfFilterStatus::MIPFilter)))
    {
        const OnboardFilterTds::ObfMipStatus* mipObj = dynamic_cast<const OnboardFilterTds::ObfMipStatus*>(tdsStatus);
       LpaMipFilter *mip = new LpaMipFilter;
       mip->initialize(masterKey,cfgKey,
            cfgId,enums::Lsf::RsdState(mipObj->getState()),
            enums::Lsf::LeakedPrescaler(mipObj->getPrescalerWord()),
            version,enums::Lsf::MIP,has, prescaleFactor);
       mip->setStatusWord(mipObj->getStatusWord());
       meta.addMip(mip);
    }

    // DFC Filter is next
    if ((tdsStatus = obf.getFilterStatus(OnboardFilterTds::ObfFilterStatus::DGNFilter)))
    {
        const OnboardFilterTds::ObfDgnStatus* dgnObj = dynamic_cast<const OnboardFilterTds::ObfDgnStatus*>(tdsStatus);
       LpaDgnFilter *dgn = new LpaDgnFilter;
       dgn->initialize(masterKey,cfgKey,
            cfgId,enums::Lsf::RsdState(dgnObj->getState()),
            enums::Lsf::LeakedPrescaler(dgnObj->getPrescalerWord()),
            version,enums::Lsf::DGN,has, prescaleFactor);
       dgn->setStatusWord(dgnObj->getStatusWord());
       meta.addDgn(dgn);
    }


  };

}
