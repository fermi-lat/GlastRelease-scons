#ifndef __TKRTRACKFITALG_H
#define __TKRTRACKFITALG_H 1
/*
#include "GaudiKernel/Algorithm.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/Track/ITkrFitTool.h"

#include "GlastSvc/Reco/IKalmanParticle.h"
*/
/** 
 * @class TkrTrackFitAlg
 *
 * @brief TkrRecon Gaudi Algorithm for controlling the fit of candidate tracks. 
 *        Gaudi Tools are used to implement a particular type of track fit, in 
 *        particular allowing a match to the output of a particular pattern 
 *        recognition algorithm. This algorithm controls their creation and use. 
 *        This algorithm depends upon input from the clustering and track finding
 *        stages of TkrRecon. Results are output to the TDS class TkrFitTrack
 * 
 * @author The Tracking Software Group
 *
 * $Header$
 */
/*
class TkrTrackFitAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrTrackFitAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrTrackFitAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    // This maintains a pointer to the particular propagator needed by the track fit
    //static IKalmanParticle* m_KalParticle;
    
private:

    /// Type of fit to perform
    std::string  m_TrackFitType;

    /// Always use Sears Craftsmen tools for the job
    ITkrFitTool* m_FitTool;
};
*/
#endif
