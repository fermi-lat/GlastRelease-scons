#ifndef ROOT_CALCLUSTER_H
#define ROOT_CALCLUSTER_H

#include "CalClusterLayerData.h"
#include "CalFitParams.h"
#include "CalMSTreeParams.h"
#include "CalMomParams.h"
#include <TObject.h>
#include <TVector3.h>
#include <map>

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
     ( const std::vector<CalClusterLayerData> &,
       const CalMSTreeParams&,
       const CalFitParams&,
       const CalMomParams &,
       const std::map <std::string, double>& classprob,
       Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits ) ;
    void init
     ( const std::vector<CalClusterLayerData> &,
       const CalMSTreeParams&,
       const CalFitParams&,
       const CalMomParams &,
       const std::map <std::string, double>& classprob,
       Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits ) ;
    virtual ~CalCluster()
     {}

    const CalMSTreeParams & getMSTreeParams() const
     { return m_mstParams;}
    const CalClusterLayerData & getLayer( int i ) const
     { return m_layers[i] ; }
    const CalFitParams & getFitParams() const
     { return m_fitParams;}
    const std::map <std::string, double>& getClassesProb() const
     {return m_classesProb;}   
    const CalMomParams & getMomParams() const
     { return m_momParams ; }
    Double_t getRmsLong() const
     { return m_momParams.getLongRms(); }
    Double_t getRmsLongAsym() const
     { return m_momParams.getLongRmsAsym(); }
    Double_t getRmsTrans() const
     { return m_momParams.getTransRms(); }
    Double_t getSkewnessLong() const
     { return m_momParams.getLongSkewness(); }
    Int_t getNumSaturatedXtals() const
     { return m_numSaturatedXtals; }
    Int_t getNumTruncXtals() const
     { return m_numTruncXtals ; }
    UInt_t getStatusBits() const
     { return m_statusBits ; }
     
    /// Access individual probabilities
    Double_t getTopologyProb(std::string) const;
    Double_t getGamProb()    const {return getTopologyProb("gam");}
    
    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalCluster &, const std::string & name = "" ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    CalClusterLayerData m_layers[ROOT_NUMCALLAYERS] ;
    CalMomParams m_momParams ;
    Int_t m_numSaturatedXtals;
    Int_t m_numTruncXtals ;
    UInt_t m_statusBits ;
    CalFitParams m_fitParams;
    CalMSTreeParams m_mstParams;    //! Parameters of the "Minimum Spanning Tree"
    std::map <std::string, double> m_classesProb;     // Classification probabilities

    ClassDef(CalCluster,8)
};


#endif 


