#ifndef RootConvert_AdfDigiConvert_H
#define RootConvert_AdfDigiConvert_H 1

#include "AncillaryDataEvent/Digi.h"

#include "digiRootData/TaggerHit.h"
#include "digiRootData/QdcHit.h"
#include "digiRootData/AdfDigi.h"


namespace RootPersistence {

 
  /// Fill persitent verison of QdcHit
  void convert( const AncillaryData::QdcHit*, QdcHit&);
  /// Fill transient verison of QdcHit
  void convert ( const QdcHit*, AncillaryData::QdcHit&);

  /// Fill persitent verison of TaggerHit
  void convert( const AncillaryData::TaggerHit*, TaggerHit&);
  /// Fill transient verison of TaggerHit
  void convert ( const TaggerHit*, AncillaryData::TaggerHit&);

  /// Fill persitent verison of AdfDigi
  void convert( const AncillaryData::Digi&, AdfDigi& ) ;
  /// Fill transient verison of AdfDigi
  void convert( const AdfDigi&, AncillaryData::Digi& ) ;

 

}

#endif












