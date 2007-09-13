#ifndef RootConvert_OnboardFilterConvert_H
#define RootConvert_OnboardFilterConvert_H 1

#include "OnboardFilterTds/FilterStatus.h"
#include "OnboardFilterTds/ObfFilterStatus.h"
#include "EFC/TFC_prjDef.h"

#include "digiRootData/FilterStatus.h"
#include "digiRootData/ObfFilterStatus.h"

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

  /// New ObfFilterStatus conversion routines begin here
  /// Fill persistent verison of ObfFilterStatus 
  void convert( const OnboardFilterTds::ObfFilterStatus&, ObfFilterStatus& ) ;
  /// Fill transient verison of ObfFilterStatus 
  void convert( const ObfFilterStatus&, OnboardFilterTds::ObfFilterStatus& ) ;

  /// Fill persistent verison of ObfGammaStatus 
  void convert( const OnboardFilterTds::IObfStatus&, ObfGammaStatus& ) ;
  /// Fill transient verison of ObfFilterStatus 
  void convert( const IObfStatus&, OnboardFilterTds::ObfGammaStatus& ) ;

  /// Fill persistent verison of ObfHFCStatus 
  void convert( const OnboardFilterTds::IObfStatus&, ObfHFCStatus& ) ;
  /// Fill transient verison of ObfFilterStatus 
  void convert( const IObfStatus&, OnboardFilterTds::ObfHFCStatus& ) ;

  /// Fill persistent verison of ObfMipStatus 
  void convert( const OnboardFilterTds::IObfStatus&, ObfMipStatus& ) ;
  /// Fill transient verison of ObfFilterStatus 
  void convert( const IObfStatus&, OnboardFilterTds::ObfMipStatus& ) ;

  /// Fill persistent verison of ObfDFCStatus 
  void convert( const OnboardFilterTds::IObfStatus&, ObfDFCStatus& ) ;
  /// Fill transient verison of ObfFilterStatus 
  void convert( const IObfStatus&, OnboardFilterTds::ObfDFCStatus& ) ;

}

#endif












