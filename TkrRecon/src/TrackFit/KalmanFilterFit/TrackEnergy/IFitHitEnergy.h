/**
 * @class IFitHitEnergy
 *
 * @brief Interface class definition for assigning energy to hits during the generic Kalman Filter Fit
 *        It is assumed that the energies stored in GLAST tracks is the Kinetic Energy of the particle
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

    /// Initialize the hit energy before doing track fit
    virtual double initialHitEnergy(const Event::TkrTrack&    patCand, 
                                    const Event::TkrTrackHit& candHit, 
                                    const double              trkEnergy) = 0;
    
    /// Updates the energy 
    virtual double updateHitEnergy(const double curEnergy, const double radLen) = 0;

    /// Retrieve the kinetic energy 
    virtual double getHitEnergy(const double energy) = 0;

    /// Convert kinetic energy to p * beta 
    virtual double kinETopBeta(const double energy) = 0;

    /// Convert p * beta to kinetic energy
    virtual double pBetaToKinE(const double energy) = 0;
};


#endif
