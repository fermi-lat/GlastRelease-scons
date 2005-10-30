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
        , EvtEnergyCorr("EvtEnergyCorr")
        , McEnergy("McEnergy")
    {}

    //function to generate goodcal test

    virtual bool isgood()
    {
        double energyRatio = EvtEnergyCorr/McEnergy;
        return fabs(energyRatio-1)<0.35;
    }

private:
    Entry EvtEnergyCorr;
    Entry McEnergy;

};

