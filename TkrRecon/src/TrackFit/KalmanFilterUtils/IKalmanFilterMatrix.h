/**
 * @class IKalmanFilterMatrix
 *
 * @brief Defines an interface for matrices used in the Kalman Filter fit
 *        This interface allows for different implementations of the main fit matrices 
 *        (Transport, Projection, Process Noise) allowing the fit to change when needed
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef IKalmanFilterMatrix_h
#define IKalmanFilterMatrix_h

#include "KalmanFilterDefs.h"

class KalmanFilterInit;
namespace Event
{
    class TkrTrackHit;
}
namespace idents {class TkrId;};

class IKalmanFilterMatrix 
{
public:
    // Define virtual methods for returning a matrix used in the Kalman Filter Fit
    virtual KFmatrix& operator()(const double &deltaZ) = 0;
    virtual KFmatrix& operator()(const idents::TkrId &id) = 0;
    virtual KFmatrix& operator()(const Event::TkrTrackHit& referenceHit,
                                 const Event::TkrTrackHit& filterHit, 
                                 const double&             eStart, 
                                 bool                      forward = true) = 0;
};


#endif
