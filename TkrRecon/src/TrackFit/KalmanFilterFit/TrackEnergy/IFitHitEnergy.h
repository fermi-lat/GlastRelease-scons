/**
 * @class IFitHitEnergy
 *
 * @brief Interface class definition for assigning energy to hits during the generic Kalman Filter Fit
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef IFitHitEnergy_h
#define IFitHitEnergy_h

namespace Event
{
    class TkrTrack;
    class TkrTrackHit;
}

class IFitHitEnergy 
{
public:

    virtual double initialHitEnergy(const Event::TkrTrack&    patCand, 
                                    const Event::TkrTrackHit& candHit, 
                                    const double              trkEnergy) = 0;
    virtual double updateHitEnergy(const double curEnergy, const double radLen) = 0;
    virtual double getHitEnergy(const double energy) = 0;
};


#endif
