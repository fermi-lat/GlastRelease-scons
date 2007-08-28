#ifndef RootConvert_OnboardFilterConvert_H
#define RootConvert_OnboardFilterConvert_H 1

#include "OnboardFilterTds/FilterStatus.h"
#include "EFC/TFC_prjDef.h"

#include "digiRootData/FilterStatus.h"

namespace RootPersistence {

  //void convert( const OnboardFilterTds::LogData*, LogData*& ) ;
  /// Fill transient verison of LogData 
  //void convert( const LogData*, OnboardFilterTds::LogData*& ) ;

  void convert(const TFC_hit &tdsObj, TfcHit& rootObj);
  void convert(const TfcHit &rootObj, TFC_hit& tdsObj);

  void convert(const TFC_prjPrms &tdsObj, TfcPrjPrms& rootObj);
  void convert(const TfcPrjPrms &rootObj, TFC_prjPrms& tdsObj);

  void convert(const TFC_prjList &tdsObj, TfcPrjList& rootObj);
  void convert(const TfcPrjList &rootObj, TFC_prjList& tdsObj);

  void convert(const TFC_prjDir &tdsObj, TfcPrjDir& rootObj);
  void convert(const TfcPrjDir &rootObj, TFC_prjDir& tdsObj);

  /// Fill persitent verison of TFC_projectionCol
  void convert( const TFC_prjs*, TfcProjectionCol&);
  /// Fill transient verison of TFC_projectionCol 
  void convert( const TfcProjectionCol&, TFC_prjs*&);

  /// Fill persitent verison of tracks 
  void convert( const std::vector<OnboardFilterTds::track>&, TObjArray&);
  /// Fill transient verison of tracks 
  void convert ( const TObjArray*, std::vector<OnboardFilterTds::track>&);

  /// Fill persitent verison of FilterStatus 
  void convert( OnboardFilterTds::FilterStatus&, FilterStatus& ) ;
  /// Fill transient verison of FilterStatus 
  void convert( const FilterStatus&, OnboardFilterTds::FilterStatus& ) ;



}

#endif












