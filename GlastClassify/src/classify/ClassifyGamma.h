/**@file ClassifyGamma.h
@brief 

*/
#pragma once
#include "GlastClassify.h"

#include <cmath>

class ClassifyGamma : public GlastClassify
{
public:
    ClassifyGamma(const std::string& info_path)
        : GlastClassify(info_path)
    {setbkgnd();}

    //functions to check or declare variables

    virtual void define(std::vector<std::string>& all_names)
    {
        m_acddocaindex = find_index("AcdDoca");
    }

    //function to generate goodcal test

    virtual bool accept()
    {
        return datum(m_acddocaindex) > 250;
    }

private:
    int m_acddocaindex;

};

