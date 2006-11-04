/**
 * @class TkrControl
 *
 * @brief A singleton class which contains control parameters 
 * for the tracking reconstruction.

 * Class is intended to be instantiated first by TkrInitSvc which 
 * has the ability to change parameters via job options parameters.
 * Once instantiated, used by pattern recognition and track fitting. 
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */


#ifndef TKRCONTROL_H
#define TKRCONTROL_H

#include <string>

//############################################
class TkrControl
//############################################
{
public:

    /// Define as a singleton object
    static TkrControl* getPtr();

    /// Retrieve values
    const double getMinEnergy()          {return m_minEnergy;         }
    const double getIniErrSlope()        {return m_iniErrorSlope;     }
    const bool   getPlaneEnergies()      {return m_planeEnergies;     }
    const bool   getTestWideClusters()   {return m_testWideClusters;  }


    /// Allow for control variables to be set at initialization
    void setMinEnergy(     double minEnergy) {m_minEnergy         = minEnergy;}
    void setIniErrSlope(    double errSlp)   {m_iniErrorSlope     = errSlp;  }
    void setPlaneEnergies(  bool   enePlane) {m_planeEnergies     = enePlane;}
    void setTestWideClusters( bool test)     {m_testWideClusters = test;}

private:
    /// private Constructor
    TkrControl();
 
    /// Pointer to the singleton object
    static TkrControl* m_this;

    /// Data members

    double m_minEnergy;          // Min tracking energy (MeV)
    double m_iniErrorSlope;      // First Hit error in Kalman: 10 deg 
    bool   m_planeEnergies;      // Decrease particle energies by exp(-rad_len)
    bool   m_testWideClusters;   // no checks on wide clusters, for heavy ions
}; 

#endif
