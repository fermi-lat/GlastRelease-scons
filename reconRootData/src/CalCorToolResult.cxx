
#include <reconRootData/CalCorToolResult.h>
#include <commonRootData/RootDataUtil.h>
#include "TParameter.h"
#include "Riostream.h"

ClassImp(CalCorToolResult)

void CalCorToolResult::init
  ( const TString & correctionName,
    UInt_t statusBits,
    Double_t chiSquare,
    const CalParams & params,
    const TObjArray & freeAttributes )
 {
  m_correctionName = correctionName ;
  m_statusBits = statusBits ;
  m_chiSquare = chiSquare ;
  m_params = params ;
  m_freeAttributes.Delete() ;
  m_freeAttributes.Clear() ;
  TIter freeParamIter(&freeAttributes) ;
  TParameter<Double_t> * freeParam ;
  while ((freeParam=(TParameter<Double_t> *)freeParamIter.Next()))
   { m_freeAttributes.Add(freeParam) ; }  
 }

void CalCorToolResult::Clear( Option_t * )
 {
  m_correctionName = "" ;
  m_statusBits = 0 ;
  m_chiSquare = 0. ;
  m_params.Clear() ;
  m_freeAttributes.Delete() ;
  m_freeAttributes.Clear() ;
 }

// dummy data, just for tests
void CalCorToolResult::Fake( Int_t ievent, UInt_t rank, Float_t randNum )
 {
  CalParams p ;
  p.Fake(ievent,rank,randNum) ;
  TObjArray freeAttributes ;
  freeAttributes.Add(new TParameter<Double_t>("Param 1",1.)) ;
  freeAttributes.Add(new TParameter<Double_t>("Param 2",2.)) ;
  init("FakeTool",3,4.,p,freeAttributes) ;
 }
 
// for tests
Bool_t CalCorToolResult::CompareInRange( const CalCorToolResult & res, const std::string & name ) const
 {
  Bool_t result = true ;
  result = rootdatautil::CompareInRange(getCorrectionName(),res.getCorrectionName(),"CorrectionName") && result ;
  result = rootdatautil::CompareInRange(getStatusBits(),res.getStatusBits(),"StatusBits") && result ;
  result = rootdatautil::CompareInRange(getChiSquare(),res.getChiSquare(),"ChiSquare") && result ;
  result = getParams().CompareInRange(res.getParams()) && result ;
  result = result &&
    rootdatautil::TObjArrayCompareInRange< TParameter<Double_t> >(
      &getFreeAttributes(),&res.getFreeAttributes(),
      "FreeAttributes"
    ) ;
//  const TObjArray & fp1 = getFreeAttributes() ;
//  const TObjArray & fp2 = res.getFreeAttributes() ;
//  result = rootdatautil::CompareInRange(fp1.GetEntries(),fp2.GetEntries(),"Number of free params")  && result ;
//  TIter fp1Iter(&fp1), fp2Iter(&fp2) ;
//  const TParameter<Double_t> * p1, * p2 ;
//  while ( ((p1=(const TParameter<Double_t> *)fp1Iter.Next())!=0) &&
//          ((p2=(const TParameter<Double_t> *)fp2Iter.Next())!=0) )
//   {
//    result = rootdatautil::CompareInRange(p1->GetName(),p2->GetName(),"Parameter Name")  && result ;
//    result = rootdatautil::CompareInRange(p1->GetVal(),p2->GetVal(),"Parameter Value")  && result ;
//   }
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;
 }

CalCorToolResult::CalCorToolResult
  ( const TString & correctionName,
    UInt_t statusBits,
    Double_t chiSquare,
    const CalParams & params,
    const TObjArray & freeAttributes )
 { init(correctionName,statusBits,chiSquare,params,freeAttributes) ; }

void CalCorToolResult::Print( Option_t * ) const
 {
  std::cout<<m_correctionName<<std::endl ;
  m_params.Print() ;
  std::cout<<"ChiSquare: "<<m_chiSquare<<std::endl ;
  m_freeAttributes.Print() ;
 }

CalCorToolResult::CalCorToolResult()
 { Clear() ; }

CalCorToolResult::~CalCorToolResult()
 { m_freeAttributes.Delete() ; }
