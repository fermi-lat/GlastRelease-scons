#ifndef RootConvert_LsfDigiConvert_H
#define RootConvert_LsfDigiConvert_H 1

#include "LdfEvent/LsfMetaEvent.h"
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
#include "digiRootData/Ccsds.h"
#include "digiRootData/LsfKeys.h"

namespace RootPersistence {

  /// Fill persitent verison of GemTime
  void convert( const lsfData::GemTime&, GemTime& ) ;
  /// Fill transient verison of GemTime
  void convert( const GemTime&, lsfData::GemTime& ) ;

  /// Fill persitent verison of TimeTone
  void convert( const lsfData::TimeTone&, TimeTone& ) ;
  /// Fill transient verison of TimeTone
  void convert( const TimeTone&, lsfData::TimeTone& ) ;

  /// Fill persitent verison of Time
  void convert( const lsfData::Time&, LsfTime& ) ;
  /// Fill transient verison of Time
  void convert( const LsfTime&, lsfData::Time& ) ;

  /// Fill persitent verison of RunInfo
  void convert( const lsfData::RunInfo&, RunInfo& ) ;
  /// Fill transient verison of RunInfo
  void convert( const RunInfo&, lsfData::RunInfo& ) ;
  
  /// Fill persitent verison of DatagramInfo
  void convert( const lsfData::DatagramInfo&, DatagramInfo& ) ;
  /// Fill transient verison of DatagramInfo
  void convert( const DatagramInfo&, lsfData::DatagramInfo& ) ;

  /// Fill persitent verison of GemScalers
  void convert( const lsfData::GemScalers&, GemScalers& ) ;
  /// Fill transient verison of GemScalers
  void convert( const GemScalers&, lsfData::GemScalers& ) ;

  /// Fill persitent verison of Configuration
  /// this handles the polymorphism of Configuration
  void convert( const lsfData::Configuration*, Configuration*& ) ;
  /// Fill transient verison of Configuration
  void convert( const Configuration*, lsfData::Configuration*& ) ;

  /// Fill persitent verison of LpaConfiguration
  void convert( const lsfData::LpaConfiguration*, LpaConfiguration&);
  /// Fill transient verison of LpaConfiguration
  void convert( const LpaConfiguration*, lsfData::LpaConfiguration*&);

  /// Fill persitent verison of LciAcdConfiguration
  void convert( const lsfData::LciAcdConfiguration*, LciAcdConfiguration&);
  /// Fill transient verison of LciAcdConfiguration
  void convert ( const LciAcdConfiguration*, lsfData::LciAcdConfiguration*&);

  /// Fill persitent verison of LciCalConfiguration
  void convert( const lsfData::LciCalConfiguration*, LciCalConfiguration&);
  /// Fill transient verison of LciCalConfiguration
  void convert ( const LciCalConfiguration*, lsfData::LciCalConfiguration*&);

  /// Fill persistent verison of LciTkrConfiguration
  void convert( const lsfData::LciTkrConfiguration*, LciTkrConfiguration&);
  /// Fill transient verison of LciTkrConfiguration
  void convert ( const LciTkrConfiguration*, lsfData::LciTkrConfiguration*&);


  /// Fill persistent verison of LpaKeys
  void convert( const lsfData::LpaKeys*, LpaKeys&);
  /// Fill transient verison of LpaKeys
  void convert ( const LpaKeys*, lsfData::LpaKeys*&);

  /// Fill persistent verison of LciKeys
  void convert( const lsfData::LciKeys*, LciKeys&);
  /// Fill transient verison of LciKeys
  void convert ( const LciKeys*, lsfData::LciKeys*&);

  /// Fill persistent verison of MetaEvent
  void convert( const LsfEvent::MetaEvent&, MetaEvent& ) ;

  /// Fill transient verison of MetaEvent
  void convert( const MetaEvent&, LsfEvent::MetaEvent& ) ;

  /// Fill persistent verison of Ccsds 
  void convert (const LsfEvent::LsfCcsds&, Ccsds&);
  /// Fill transient verison of Ccsds 
  void convert(const Ccsds&, LsfEvent::LsfCcsds&);


  /// Special case to fill Flight OBF status words when running MC
  void convert( const OnboardFilterTds::ObfFilterStatus&, MetaEvent& );
}

#endif












