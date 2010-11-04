#ifndef ROOT_CALCLUSTER_H
#define ROOT_CALCLUSTER_H

#include "CalClusterLayerData.h"
#include "CalFitParams.h"
#include "CalMSTreeParams.h"
#include "CalParams.h"
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
       const CalParams &,
       const std::map <std::string, double>& classprob,          
       Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans, Double_t skewLong,
       Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits ) ;
    void init
     ( const std::vector<CalClusterLayerData> &,
       const CalMSTreeParams&,
       const CalFitParams&,
       const CalParams &,
       const std::map <std::string, double>& classprob,                
       Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans, Double_t skewLong,
       Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits ) ;
    virtual ~CalCluster()
     {}

    const CalMSTreeParams & getMSTreeParams() const
     { return m_mstreeParams;}
    const CalClusterLayerData & getLayer( int i ) const
     { return m_layers[i] ; }
    const CalFitParams & getFitParams() const
     { return m_fitParams;}
    const std::map <std::string, double>& getClassesProb() const
     {return m_classesProb;}   
    const CalParams & getParams() const
     { return m_params ; }     
    Double_t getRmsLong() const
     { return m_rmsLong ; }
    Double_t getRmsLongAsym() const
     { return m_rmsLongAsym ; }
    Double_t getRmsTrans() const
     { return m_rmsTrans ; }
    Double_t getSkewnessLong() const
     { return m_skewnessLong ; }
    Int_t getNumSaturatedXtals() const
     { return m_numSaturatedXtals;}
    Int_t getNumTruncXtals() const
     { return m_numTruncXtals ; }
    UInt_t  getStatusBits() const
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
    CalParams m_params ;
    Double32_t m_rmsLong ;
    Double32_t m_rmsLongAsym ;
    Double32_t m_rmsTrans ;
    Double32_t m_skewnessLong ;
    Int_t m_numSaturatedXtals;
    Int_t m_numTruncXtals ;
    UInt_t m_statusBits ;
    CalFitParams m_fitParams;
    CalMSTreeParams m_mstreeParams;    //! Parameters of the "Minimum Spanning Tree"
    std::map <std::string, double> m_classesProb;     // Classification probabilities

    ClassDef(CalCluster,7)
};


#endif 


