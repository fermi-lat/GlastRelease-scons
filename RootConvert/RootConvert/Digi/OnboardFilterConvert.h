#ifndef RootConvert_OnboardFilterConvert_H
#define RootConvert_OnboardFilterConvert_H 1

#include "OnboardFilter/FilterStatus.h"
#include "EFC/TFC_projectionDef.h"

#include "digiRootData/FilterStatus.h"

namespace RootPersistence {

  //void convert( const OnboardFilterTds::LogData*, LogData*& ) ;
  /// Fill transient verison of LogData 
  //void convert( const LogData*, OnboardFilterTds::LogData*& ) ;

  void convert(const TFC_projectionDir &tdsObj, TfcProjectionDir& rootObj);
  void convert(const TfcProjectionDir &rootObj, TFC_projectionDir& tdsObj);

  /// Fill persitent verison of TFC_projectionCol
  void convert( const TFC_projections*, TfcProjectionCol&);
  /// Fill transient verison of TFC_projectionCol 
  void convert( const TfcProjectionCol&, TFC_projections*&);

  /// Fill persitent verison of tracks 
  void convert( const std::vector<OnboardFilterTds::track>&, TObjArray&);
  /// Fill transient verison of tracks 
  void convert ( const TObjArray*, std::vector<OnboardFilterTds::track>&);

  /// Fill persitent verison of FilterStatus 
  void convert( OnboardFilterTds::FilterStatus&, FilterStatus& ) ;
  /// Fill transient verison of FilterStatus 
  void convert( const FilterStatus&, OnboardFilterTds::FilterStatus*& ) ;



}

#endif












