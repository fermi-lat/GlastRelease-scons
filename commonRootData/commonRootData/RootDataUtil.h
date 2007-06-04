
#ifndef ROOTDATAUTIL_H
#define ROOTDATAUTIL_H

#include <commonRootData/idents/AcdId.h>
#include <commonRootData/idents/CalXtalId.h>
#include <commonRootData/idents/TkrId.h>
#include <commonRootData/idents/TowerId.h>
#include <commonRootData/idents/VolumeIdentifier.h>
#include <commonRootData/Exception.h>

#include "Rtypes.h"
#include "TObjArray.h"
#include "TParameter.h"
class TVector3 ;
class TLorentzVector ;
class TString ;

#include <string>
#include <vector>

namespace rootdatautil {
    
    extern void CheckBounds(long nIndex, long nLower, long nUpper);
    extern void CheckBufSize(unsigned long nBufSize, unsigned long nCorrectSize); 

    // low level comparisons
    bool CompareInRange( Float_t, Float_t, const std::string & name ="unknown Double_t" ) ;
    bool CompareInRange( Double_t, Double_t, const std::string & name ="unknown Double_t" ) ;
    bool CompareInRange( Short_t, Short_t, const std::string & name ="unknown integer" ) ;
    bool CompareInRange( Int_t, Int_t, const std::string & name ="unknown integer" ) ;
    bool CompareInRange( UInt_t, UInt_t, const std::string & name ="unknown integer" ) ;
    bool CompareInRange( ULong_t, ULong_t, const std::string & name ="unknown integer" ) ;
    bool CompareInRange( Bool_t, Bool_t, const std::string & name ="unknown bool" ) ;
    bool CompareInRange( Char_t, Char_t, const std::string & name ="unknown Char_t" ) ;
    bool CompareInRange( const char *, const char *, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const TString &, const TString &, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const std::string &, const std::string &, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const TVector3 &, const TVector3 &, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const TLorentzVector &, const TLorentzVector &, const std::string & name ="unknown lorentz vector" ) ;
    bool CompareInRange( void *, void *, const std::string & name ="unknown pointer" ) ;
    bool CompareInRange( const TParameter<Double_t> &, const TParameter<Double_t> &, const std::string & name ="unknown TParameter<Double_t>" ) ;
    bool CompareInRange( const AcdId &, const AcdId &, const std::string & name ="unknown AcdId" ) ;
    bool CompareInRange( const CalXtalId &, const CalXtalId &, const std::string & name ="unknown CalXtalId" ) ;
    bool CompareInRange( const CalXtalId::CalTrigMode &, const CalXtalId::CalTrigMode &, const std::string & name ="unknown CalXtalId::CalTrigMode" ) ;
    bool CompareInRange( const commonRootData::TkrId &, const commonRootData::TkrId &, const std::string & name ="unknown TkrId" ) ;
    bool CompareInRange( const TowerId &, const TowerId &, const std::string & name ="unknown TowerId" ) ;
    bool CompareInRange( const VolumeIdentifier &, const VolumeIdentifier &, const std::string & name ="unknown VolumeIdentifier" ) ;

    // adapters
    template <class T>
    bool CompareInRange( const T & obj1, const T & obj2, const std::string & name = "" )
     {
      if (name=="") return obj1.CompareInRange(obj2) ;
      else  return obj1.CompareInRange(obj2,name) ;
     }
    template <class T>
    bool CompareInRange( const T * obj1, const T * obj2, const std::string & name = "" )
     { return CompareInRange((void*)obj1,(void*)obj2,name) ; }

    // std collections
    bool CompareInRange(
      const std::vector<Double_t> &,
      const std::vector<Double_t> &,
      const std::string & name ="unknown std::vector<Double_t>" ) ;
    bool CompareInRange(
      const std::vector<UInt_t> &,
      const std::vector<UInt_t> &,
      const std::string & name ="unknown std::vector<UInt_t>" ) ;
    bool CompareInRange(
      const std::vector<AcdId> &,
      const std::vector<AcdId> &,
      const std::string & name ="unknown std::vector<AcdId>" ) ;

    // root collections
    template <class Element>
    bool TObjArrayCompareInRange(
        const TObjArray * col1,
        const TObjArray * col2,
        const std::string & name = "unknown TObjArray" ) {
        
        bool res = CompareInRange(col1->GetEntries(),col2->GetEntries(),name+" Size") ;
        if (!res) return res ;
         
        TIter rootIt1(col1);
        TIter rootIt2(col2);
        Element * element1 ;
        Element * element2 ;
        while ( (element1 = (Element*)rootIt1.Next()) &&
                (element2 = (Element*)rootIt2.Next()) ) {
            res = res && CompareInRange(*element1,*element2,name+" Element") ;
        }

        return res ;
    }
} ;


#endif
