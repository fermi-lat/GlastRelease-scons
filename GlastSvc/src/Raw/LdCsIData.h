// LdCsIData.h: interface for the LdCsIData class.
//
//////////////////////////////////////////////////////////////////////
#ifndef LDCSI_DATA_H
#define LDCSI_DATA_H 1



#include "GlastEvent/data/TdCsIData.h"
#include "instrument/CsIDetector.h"


/*! Takes care of loading a TdCsIData object from CsIDetector. This 
    class was created to remove the dependence of GlastEvent on 
    instrument.

    Author: Ian Gable
*/
class LdCsIData : public TdCsIData {

public:
    LdCsIData(int numLay); 
    ~LdCsIData(){};
    
    //!load local array of xtal objects from a CsIDetector
    void load (const CsIDetector& xtal, idents::ModuleId tower);
    
    
};

#endif