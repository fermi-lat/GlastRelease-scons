#ifndef RootConvert_TkrTruncationInfoConvert_H
#define RootConvert_TkrTruncationInfoConvert_H 1

//needed because CLID and 
//InterfaceID are not prototyped in TkrTruncatedData
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/DataObject.h"

#include <Event/Recon/TkrRecon/TkrTruncationInfo.h>
#include <reconRootData/TkrTruncationData.h>

//#include <pair>
#include <utility>
#include <TCollection.h>

namespace RootPersistence {

void convert(const std::pair<Event::SortId, Event::TkrTruncatedPlane> &, TkrTruncationData &);
void convert( const Event::TkrTruncationInfo &, TkrTruncationData & ) ;
void convert( const TCollection &, Event::TkrTruncationInfo & ) ;

}

#endif
