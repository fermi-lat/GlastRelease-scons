/**
 * @class McLayerHit
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "TkrRecon/MonteCarlo/McLayerHit.h"

namespace Event {
McLayerHit::McLayerHit(const Event::McParticle* particle): 
m_McParticle(particle), m_cluster(0)
{
    m_PositionHitsVec.clear();
}

McLayerHit::~McLayerHit()
{
    return;
}

void McLayerHit::addMcPositionHit(const Event::McPositionHit* posHit)
{
    if (m_PositionHitsVec.size() == 0)
    {
        m_volIdent = posHit->volumeID();
    }

    m_PositionHitsVec.push_back(posHit);

    return;
}
};
