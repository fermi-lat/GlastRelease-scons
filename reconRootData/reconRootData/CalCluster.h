#ifndef ROOT_CALCLUSTER_H
#define ROOT_CALCLUSTER_H

#include "CalClusterLayerData.h"
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
     ( const std::vector<CalClusterLayerData> &, const CalParams &,
       Double_t rmsLong, Double_t rmsTrans ) ;
    void init
     ( const std::vector<CalClusterLayerData> &, const CalParams &,
       Double_t rmsLong, Double_t rmsTrans ) ;
    virtual ~CalCluster()
     {}

    const CalClusterLayerData & getLayer( int i ) const
     { return m_layers[i] ; }
    const CalParams & getParams() const
     { return m_params ; }
    Double_t getRmsLong() const
     { return m_rmsLong ; }
    Double_t getRmsTrans() const
     { return m_rmsTrans ; }

    void Clear( Option_t * option ="" ) ;
    void Print( Option_t * option ="" ) const ;

private:

    CalClusterLayerData m_layers[ROOT_NUMCALLAYERS] ;
    CalParams m_params ;
    Double_t m_rmsLong;
    Double_t m_rmsTrans;

    ClassDef(CalCluster,5)
};


#endif 


