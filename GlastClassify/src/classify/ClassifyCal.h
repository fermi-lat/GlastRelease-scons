/**@file ClassifyCal.h
@brief 

*/
#include "GlastClassify.h"

#include <cmath>

class ClassifyCal : public GlastClassify
{
public:
    ClassifyCal(const std::string& info_path)
        : GlastClassify(info_path)
    {}

    //functions to check or declare variables

    virtual void define(std::vector<std::string>& all_names)
    {
        m_energyIndex   =find_index("EvtEnergyCorr");
        m_calEnergyIndex= find_index("CalEnergyRaw");
        m_CalTotRLnIndex = find_index("CalTotRLn");
        m_mcEnergy = subdefine(all_names, "McEnergy");
    }

    //function to generate goodcal test

    virtual bool isgood()
    {
        double energyRatio = datum(m_energyIndex)/datum(m_mcEnergy);
        return fabs(energyRatio-1)<0.35;
    }

    virtual bool accept()
    {
        return  datum(m_calEnergyIndex) > 5.0 && datum(m_CalTotRLnIndex)>4.0;

    }
    float calEnergy()const{ return datum(m_calEnergyIndex); }

private:
    int m_energyIndex;
    int m_mcEnergy;
    int m_calEnergyIndex;
    int m_CalTotRLnIndex;
};

