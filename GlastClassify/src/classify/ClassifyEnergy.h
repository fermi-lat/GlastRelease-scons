/**@file ClassifyEnergy.h
@brief Define ClassifyEnergy

$Header$
*/
#include "GlastClassify.h"

#include <cmath>

/** @class GlastClassify
    @brief Manage classification of the various energy estimates

*/

class ClassifyEnergy : public GlastClassify
{
public:
    typedef enum {PARAM, LASTLAYER, PROFILE, TRACKER, BAD} Type;
    Type m_type;

    ClassifyEnergy(const std::string& info_path)
        : GlastClassify(info_path)
        , EvtEnergyCorr("EvtEnergyCorr")
        , McEnergy("McEnergy")
        , McLogEnergy("McLogEnergy")
        , CalEnergyRaw("CalEnergyRaw")
        , CalTklEnergy("CalTklEnergy")
        , CalCfpEnergy("CalCfpEnergy")
        , CalLllEnergy("CalLllEnergy")
        , CalTotRLn("CalTotRLn")
    { 
        if( info_path.find("param")!=std::string::npos) m_type=PARAM;
        else if( info_path.find("lastlayer")!=std::string::npos) m_type=LASTLAYER;
        else if( info_path.find("profile")!=std::string::npos) m_type=PROFILE;
        else if( info_path.find("tracker")!=std::string::npos) m_type=TRACKER;
        else {
            m_type = BAD;
        }
    }

    //function to generate good eneregy measurement test

    virtual bool isgood()
    {
        double energyRatio = m_energy/McEnergy;
        double tolerance = 2.* energyResModel();
        return fabs(energyRatio-1) < tolerance; 
    }

    virtual bool accept()
    {
        if(CalEnergyRaw<5. || CalTotRLn<4.) return false; 	 
 
        switch( m_type )
        {
        case PARAM: m_energy = EvtEnergyCorr; break;
        case PROFILE: m_energy = CalCfpEnergy; break;
        case TRACKER: m_energy = CalTklEnergy; break;
        case LASTLAYER: m_energy = CalLllEnergy; break;
        default: m_energy=0;
        };
        return m_energy>0.;
    }

    double energyResModel(){
        // defined by Atwood 
        return 0.05 + 0.72/pow(McLogEnergy,3);
    }
private:
    double m_energy;
    Entry EvtEnergyCorr;
    Entry McEnergy;
    Entry McLogEnergy;
    Entry CalEnergyRaw;
    Entry CalTklEnergy;
    Entry CalCfpEnergy;
    Entry CalLllEnergy;
    Entry CalTotRLn;
};

