#ifndef ROOT_CALCLUSTER_H
#define ROOT_CALCLUSTER_H

#include "CalClusterLayerData.h"
#include "CalFitParams.h"
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

#define ROOT_NUMCALLAYERS 8

/** @class CalCluster
 * @brief Root object that stores a cluster of CsI detections
 *
 * Provides access to the positions and energy depositions for a
 * cluster of CsI logs 
 *
 * @author A. Chekhtman
 *
 * $Header$
 */
class CalCluster : public TObject
{ 
public:
    
    CalCluster()
     { Clear() ; }
    CalCluster
     ( const std::vector<CalClusterLayerData> &, const CalFitParams&, const CalParams &,
       Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans,
       Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits ) ;
    void init
     ( const std::vector<CalClusterLayerData> &, const CalFitParams&, const CalParams &,
       Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans,
       Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits ) ;
    virtual ~CalCluster()
     {}

    const CalClusterLayerData & getLayer( int i ) const
     { return m_layers[i] ; }
    const CalFitParams & getFitParams() const
     { return m_fitParams;}
    const CalParams & getParams() const
     { return m_params ; }
    Double_t getRmsLong() const
     { return m_rmsLong ; }
    Double_t getRmsLongAsym() const
     { return m_rmsLongAsym ; }
    Double_t getRmsTrans() const
     { return m_rmsTrans ; }
    Int_t getNumSaturatedXtals() const
     { return m_numSaturatedXtals;}
    Int_t getNumTruncXtals() const
     { return m_numTruncXtals ; }
    UInt_t  getStatusBits() const
     { return m_statusBits ; }
     
    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalCluster &, const std::string & name = "" ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    CalClusterLayerData m_layers[ROOT_NUMCALLAYERS] ;
    CalFitParams m_fitParams;
    CalParams m_params ;
    Double32_t m_rmsLong ;
    Double32_t m_rmsLongAsym ;
    Double32_t m_rmsTrans ;
    Int_t m_numSaturatedXtals;
    Int_t m_numTruncXtals ;
    UInt_t m_statusBits ;

    ClassDef(CalCluster,6)
};


#endif 


