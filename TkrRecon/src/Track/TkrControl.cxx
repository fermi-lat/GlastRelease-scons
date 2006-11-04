// File and Version Information:
//      $Header$
//
// Description:
//      Implements singleton class for storing and retrieving 
//      tracker recon control parameters.
//
// Author:
//      The Tracking Software Group  


#include "src/Track/TkrControl.h"

TkrControl* TkrControl::m_this = 0;

TkrControl::TkrControl()
{
    m_minEnergy          = 30.0; // Min tracking energy (MeV)
    m_iniErrorSlope      = 0.17; // First Hit error in Kalman: 10 deg 
    m_planeEnergies      = true; // Decrease particle energies by exp(-rad_len)
    m_testWideClusters   = true; // turn off for heavy ions

    return;
}

TkrControl* TkrControl::getPtr()
{
    if (!m_this) m_this = new TkrControl();

    return m_this;
}
