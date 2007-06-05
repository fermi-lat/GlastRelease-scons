#include <LdfEvent/LsfMetaEvent.h>
#include "LdfEvent/LsfCcsds.h"

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
    rootObj.initialize(tdsObj.platform(),tdsObj.dataOrigin(),tdsObj.id(),tdsObj.startTime());
  };

  void convert( const RunInfo& rootObj, lsfData::RunInfo& tdsObj) {
    tdsObj.set(rootObj.platform(),rootObj.dataOrigin(),rootObj.id(),rootObj.startTime());
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
          tdsObj->biasDac(), tdsObj->holdDelay(), triggerRoot, channelRoot);

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
         rootObj->biasDac(), rootObj->holdDelay(), triggerTds, chTds);
      tdsObj->setSoftwareKey(rootObj->softwareKey());
      tdsObj->setWriteCfg(rootObj->writeCfg());
      tdsObj->setReadCfg(rootObj->readCfg());
      tdsObj->setPeriod(rootObj->period());
      tdsObj->setFlags(rootObj->flags());
  };

  void convert( const lsfData::LciCalConfiguration* tdsObj,
                LciCalConfiguration& rootObj) {
      lsfData::LciCalConfiguration::CalTrigger trigTds = tdsObj->trigger();

      LciCalConfiguration::CalTrigger trigRoot(trigTds.le(), trigTds.he());

      lsfData::Channel chTds = tdsObj->channel();
      Channel chRoot(chTds.single(), chTds.all(), chTds.latc());

      rootObj.initialize(tdsObj->uld(), tdsObj->injected(), tdsObj->delay(),
              tdsObj->threshold(),trigRoot,chRoot);
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
          trigRoot.he());
      Channel chRoot = rootObj->channel();
      lsfData::Channel chTds(chRoot.single(), chRoot.all(), chRoot.latc());
      tdsObj->set(rootObj->uld(), rootObj->injected(), rootObj->delay(),
              rootObj->threshold(),trigTds,chTds);
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
              tdsObj->threshold(), chRoot);
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
              rootObj->threshold(), chTds);
      tdsObj->setSoftwareKey(rootObj->softwareKey());
      tdsObj->setWriteCfg(rootObj->writeCfg());
      tdsObj->setReadCfg(rootObj->readCfg());
      tdsObj->setPeriod(rootObj->period());
      tdsObj->setFlags(rootObj->flags());
  };

  void convert( const LpaKeys* rootObj,
                lsfData::LpaKeys*& tdsObj){

    tdsObj->setCDM_keys(rootObj->CDM_keys());
    tdsObj->setLATC_master(rootObj->LATC_master());
    tdsObj->setLATC_ignore(rootObj->LATC_ignore());

  };

  void convert( const lsfData::LpaKeys* tdsObj,
                LpaKeys& rootObj) {

    rootObj.setCDM_keys(tdsObj->CDM_keys());
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
  };

  void convert( const LsfEvent::LsfCcsds& tdsObj, Ccsds& rootObj) {
      rootObj.initialize(tdsObj.scid(), tdsObj.apid(), tdsObj.utc());
  };

  void convert (const Ccsds& rootObj, LsfEvent::LsfCcsds& tdsObj) {
      tdsObj.initialize(rootObj.getScid(), rootObj.getApid(), rootObj.getUtc());
  };

}
