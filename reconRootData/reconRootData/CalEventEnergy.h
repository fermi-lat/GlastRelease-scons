#ifndef ROOT_CalEventEnergy_H
#define ROOT_CalEventEnergy_H

#include "CalCorToolResult.h"
#include "CalParams.h"
#include <TObject.h>
#include <TVector3.h>

#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif


class CalEventEnergy : public TObject
{ 
public:
//    
//    CalEventEnergy()
//     { Clear() ; }
//    CalEventEnergy
//     ( const std::vector<CalClusterLayerData> &, const CalParams &,
//       Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans,
//       Int_t numTruncXtals, UInt_t statusBits ) ;
//    void init
//     ( const std::vector<CalClusterLayerData> &, const CalParams &,
//       Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans,
//       Int_t numTruncXtals, UInt_t statusBits ) ;
//    virtual ~CalEventEnergy()
//     {}
//
//    const CalClusterLayerData & getLayer( int i ) const
//     { return m_layers[i] ; }
    const CalParams & getParams() const
     { return m_params ; }
//    Double_t getRmsLong() const
//     { return m_rmsLong ; }
//    Double_t getRmsLongAsym() const
//     { return m_rmsLongAsym ; }
//    Double_t getRmsTrans() const
//     { return m_rmsTrans ; }
//    Int_t getNumTruncXtals() const
//     { return m_numTruncXtals ; }
//    UInt_t  getStatusBits() const
//     { return m_statusBits ; }
//     
    void Clear( Option_t * option ="" ) ;
    void Fake( UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalEventEnergy & ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

//    CalClusterLayerData m_layers[ROOT_NUMCALLAYERS] ;
    CalParams m_params ;
//    Double_t m_rmsLong ;
//    Double_t m_rmsLongAsym ;
//    Double_t m_rmsTrans ;
//    Int_t m_numTruncXtals ;
//    UInt_t m_statusBits ;

    ClassDef(CalEventEnergy,1)
};


#endif 


