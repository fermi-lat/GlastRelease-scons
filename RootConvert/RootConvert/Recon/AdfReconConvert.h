#ifndef RootConvert_AdfReconConvert_H
#define RootConvert_AdfReconConvert_H 1

#include "AncillaryDataEvent/Recon.h"
#include "reconRootData/AdfRecon.h"


namespace RootPersistence {

/* 
  /// Fill persitent verison of TaggerHit
  void convert( const AncillaryData::TaggerHit*, commonRootData::TaggerHit&);
  /// Fill transient verison of TaggerHit
  void convert ( const commonRootData::TaggerHit*, AncillaryData::TaggerHit&);
*/

  /// Fill persitent verison of TaggeCluster
  void convert( const AncillaryData::TaggerCluster, reconRootData::TaggerCluster*&);
  /// Fill transient verison of TaggerHit
  void convert ( const reconRootData::TaggerCluster*, AncillaryData::TaggerCluster&);

  /// Fill persitent verison of AdfRecon
  void convert( const AncillaryData::Recon&, reconRootData::AdfRecon& ) ;
  /// Fill transient verison of AdfRecon
  void convert( const reconRootData::AdfRecon&, AncillaryData::Recon& ) ;

 

}

#endif












