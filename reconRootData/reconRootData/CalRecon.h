#ifndef ROOT_CALRECON_H
#define ROOT_CALRECON_H

#include "TObject.h"
#include "TObjArray.h"

//#include <vector>
//#ifdef WIN32
//using namespace std;
//#endif

#include "CalXtalRecData.h"
#include "CalCluster.h"
#include "CalMipTrack.h"

/** @class CalRecon
 * @brief Primary Root object containing CAL reconstruction data.
 *
 * Portal through which one gains access to the CAL recon data. 
 * Currently there are 2 lists available:
 * - A collection of reconstructed CAL crystals (CalXtalRecData)
 * - A collection of CAL clusters (CalCluster)
 *
 * $Header$
*/

class CalRecon : public TObject
{

public:

    CalRecon();

    virtual ~CalRecon();

    void initialize();

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    /// provide access to list of CalClusters
    TObjArray* getCalClusterCol()  { return m_clusterCol; };

    /// provide access to the list of CalXtalRecData
    TObjArray* getCalXtalRecCol()  { return m_xtalRecCol; };

    /// provide access to the list of CalMipTracks
    TObjArray* getCalMipTrackCol() { return m_mipTrackCol; };

    void addCalCluster(CalCluster* cluster)   { m_clusterCol->Add(cluster); };
    void addXtalRecData(CalXtalRecData *xtal) { m_xtalRecCol->Add(xtal); };
    void addCalMipTrack(CalMipTrack* track)   { m_mipTrackCol->Add(track); };

private:
    /// list of CalXtalRecData
    TObjArray* m_xtalRecCol;

    /// list of CalClusters
    TObjArray* m_clusterCol;

    /// list of CalMipTracks
    TObjArray* m_mipTrackCol;
  
    ClassDef(CalRecon,3)
};

#endif
