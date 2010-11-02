// File and Version information:
// $Header:
//
//  Implementation file of Rot clone for CalMSTreeParams
//  
// Authors:
//
//        Luca Baldini, Johan Bregeon, Carmelo Sgro'
//
//

#include <reconRootData/CalMSTreeParams.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalMSTreeParams)

CalMSTreeParams::CalMSTreeParams
( double totalEnergy,
  double maxXtalEnergy,  int	numEdges,	    					  
  double minEdgeLength,  double maxEdgeLength,	    					  
  double meanEdgeLength, double meanEdgeLengthTrunc,					  
  double rmsEdgeLength,  double rmsEdgeLengthTrunc)
{
    init
    ( totalEnergy, maxXtalEnergy, numEdges, minEdgeLength, maxEdgeLength,
      meanEdgeLength, meanEdgeLengthTrunc, rmsEdgeLength, rmsEdgeLengthTrunc
    );
}

// Note -- Clear instead of clear in original class
void CalMSTreeParams::Clear(  Option_t * )
{
    init
    ( 0.,0.,0 ,0.,0.,
      0.,0.,0.,0. );
}


void CalMSTreeParams::init
( Double_t totalEnergy,
  Double_t maxXtalEnergy,  Int_t    numberOfEdges,			
  Double_t minEdgeLength,  Double_t maxEdgeLength,		      
  Double_t meanEdgeLength, Double_t meanEdgeLengthTrunc,	      
  Double_t rmsEdgeLength,  Double_t rmsEdgeLengthTrunc)
{
    m_totalEnergy    = totalEnergy;
    m_maxXtalEnergy  = maxXtalEnergy;
    m_numberOfEdges  = numberOfEdges;
    m_minEdgeLength  = minEdgeLength;
    m_maxEdgeLength  = maxEdgeLength;
    m_meanEdgeLength = meanEdgeLength;
    m_rmsEdgeLength  = rmsEdgeLength;
    m_meanEdgeLengthTrunc = meanEdgeLengthTrunc;
    m_rmsEdgeLengthTrunc  = rmsEdgeLengthTrunc;
}

// Dump the content of the container: the MST tree parameters
void CalMSTreeParams::Print( Option_t * ) const
{
    std::cout
      << m_totalEnergy	       << " " << m_maxXtalEnergy << " " << m_numberOfEdges << " "
      << m_minEdgeLength       << " " << m_maxEdgeLength << " " << m_meanEdgeLength<< " "
      << m_meanEdgeLengthTrunc << " " << m_rmsEdgeLength << " " << m_rmsEdgeLengthTrunc <<
    std::endl;
}

// create dummy data, just for tests
void CalMSTreeParams::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ )
 {
    init
    ( 1.,2.,3 ,4.,5.,
      6.,7.,8.,9. );   
 }

// CalMSTreeParams::CompareInRange, just for tests
Bool_t CalMSTreeParams::CompareInRange( const CalMSTreeParams & tp, const std::string & name ) const {
    
    Bool_t result = true ;
    
    result = rootdatautil::CompareInRange(getTotalEnergy()	  ,tp.getTotalEnergy()        ,"Total Energy") && result ;
    result = rootdatautil::CompareInRange(getMaxXtalEnergy()      ,tp.getMaxXtalEnergy()      ,"Max Xtal Energy") && result ;
    result = rootdatautil::CompareInRange(getNumberOfEdges()      ,tp.getNumberOfEdges()      ,"Number of Edges") && result ;
    result = rootdatautil::CompareInRange(getMinEdgeLength()	  ,tp.getMinEdgeLength()      ,"Min Edge Length") && result ;
    result = rootdatautil::CompareInRange(getMaxEdgeLength()	  ,tp.getMaxEdgeLength()      ,"Max Edge Length") && result ;
    result = rootdatautil::CompareInRange(getMeanEdgeLength()	  ,tp.getMeanEdgeLength()     ,"Mean of Edge Lengths") && result ;
    result = rootdatautil::CompareInRange(getMeanEdgeLengthTrunc(),tp.getMeanEdgeLengthTrunc(),"Truncated Mean of Edge Lengths") && result ;
    result = rootdatautil::CompareInRange(getRmsEdgeLength()	  ,tp.getRmsEdgeLength()      ,"Rms of Edge Lengths") && result ;
    result = rootdatautil::CompareInRange(getRmsEdgeLengthTrunc() ,tp.getRmsEdgeLengthTrunc() ,"Truncated Rms of Edge Lengths") && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for CalMSTreeParams"<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;
        
}
