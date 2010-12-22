#ifndef ROOT_CALCLUSTER_H
#define ROOT_CALCLUSTER_H

#include "CalClusterLayerData.h"
#include "CalXtalsParams.h"
#include "CalFitParams.h"
#include "CalMSTreeParams.h"
#include "CalMomParams.h"
#include "CalClassParams.h"
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
       const CalXtalsParams&,
       const CalMSTreeParams&,
       const CalFitParams&,
       const CalMomParams &,
       const CalClassParams&,
       UInt_t statusBits ) ;
    void init
     ( const std::vector<CalClusterLayerData> &,
       const CalXtalsParams&,
       const CalMSTreeParams&,
       const CalFitParams&,
       const CalMomParams&,
       const CalClassParams&,
       UInt_t statusBits ) ;
    virtual ~CalCluster()
     {}

    const CalClusterLayerData & getLayer( int i ) const
     { return m_layers[i] ; }

    std::string getProducerName() const
      { return m_producerName ; } 
    UInt_t getStatusBits() const
     { return m_statusBits ; }
    const CalXtalsParams & getXtalsParams() const
     { return m_xtalsParams;}
    const CalMSTreeParams & getMSTreeParams() const
     { return m_mstParams;}
    const CalFitParams & getFitParams() const
     { return m_fitParams;} 
    const CalMomParams & getMomParams() const
     { return m_momParams ; }
    const CalClassParams & getClassParams() const
     { return m_classParams ; }

    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalCluster &, const std::string & name = "" ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    CalClusterLayerData m_layers[ROOT_NUMCALLAYERS] ;

    std::string m_producerName ;
    UInt_t m_statusBits ;
    CalXtalsParams m_xtalsParams ;
    CalMSTreeParams m_mstParams ;
    CalFitParams m_fitParams ;
    CalMomParams m_momParams ;
    CalClassParams m_classParams ;

    ClassDef(CalCluster,9)
};


#endif 


