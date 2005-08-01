/**@file ClassifyGamma.h
@brief 

$Header$

*/
#include "GlastClassify.h"

#include <cmath>

class ClassifyGamma : public GlastClassify
{
public:
    ClassifyGamma(const std::string& info_path)
        : GlastClassify(info_path, false) // flag that not mixed.
    {
        m_calEnergyIndex = find_index("CalEnergyRaw");
        m_CalTotRLnIndex = find_index("CalTotRLn");

        m_FilterAlgStatus=  add_index( "FilterAlgStatus");
        m_acddocaindex   =  add_index( "AcdDoca");
    }

    //acceptance function, applied to background and signal

    virtual bool accept()
    {
        return datum(m_acddocaindex)   > 250
            && datum(m_calEnergyIndex) > 5.0 
            && datum(m_CalTotRLnIndex) > 4.0
            && datum(m_FilterAlgStatus)==0;
    }

private:
    int m_acddocaindex;
    int m_calEnergyIndex;
    int m_CalTotRLnIndex;
    int m_FilterAlgStatus;
};

