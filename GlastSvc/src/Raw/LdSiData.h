// LdSiData.h: interface for the LdSiData class.
//
//////////////////////////////////////////////////////////////////////
#ifndef LDSIDATA_H
#define LDSIDATA_H 1

#include "GlastEvent/data/TdSiData.h"
#include "instrument/SiDetector.h"


/*! Takes care of loading a TdSiData object from CsIDetector. This 
    class was created to remove the dependence of GlastEvent on 
    instrument.

    Author:Ian Gable
*/

class LdSiData :public TdSiData  
{
    
public:
    LdSiData(int numLay):TdSiData(numLay){
        ;
    };
    ~LdSiData(){};
    void load (const SiDetector& strips, idents::ModuleId module);
    
};

#endif 
