// $Header$
//
#ifndef MULTIPSF_H
#define MULTIPSF_H

#include "Analyze.h"
#include <vector>

class PSFanalysis;
// Analyzed multiple bins in generated energy for PSF analysis

class MultiPSF : public Analyze , public std::vector<PSFanalysis*>{
public:
    MultiPSF(const Tuple& t, char code);

    void report(std::ostream& out);

    //Put these in for use with root
    int          getListSize()        {return size();}
    PSFanalysis* getListItem(int idx) {return (*this)[idx];}
private:
    bool apply();
    double   m_bin_size; 
    std::vector<float>m_costheta_bin;
};

#endif
