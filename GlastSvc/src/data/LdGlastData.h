// LdGlastData.h: interface for the LdGlastData class.
//
//////////////////////////////////////////////////////////////////////
#ifndef LDGLAST_DATA_H
#define LDGLAST_DATA_H 1

#include "data/SiData.h"
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include "geometry/Point.h"
#include "idents/ModuleId.h"

#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"


// forward declarations
#include "data/SiData.h"
#include "data/CsIData.h"
#include "data/IVetoData.h"

#include "GlastEvent/data/TdGlastData.h"

#include <iostream>
#include <vector>

#include "instrument/MCTruth.h"

extern const CLID& CLID_LdGlastData;

//! Quick wrapper for the TdSiData and TdCsIData objects
/*! \class TdGlastData
    The TDS version of GlastData set up as a wrapped for the TDS objects.
    Only really needed as a quick fix for the tracker recon. The 
    Implementation of TdVetoData is still left, therefore take note
    of the fact that this only returns a null pointer for the Veto
    data object.

    Note that there is no path for this object in the event model and
    no converter available.
*/


class LdGlastData : virtual public TdGlastData {

public:
    LdGlastData ( const CsIData* csi, const SiData* si,const IVetoData* vd)
        : m_csiData(csi), m_siData( si), m_vetoData (vd), m_truth(MCTruth::instance()){};

    virtual const CLID& clID() const   { return LdGlastData::classID(); }
    static const CLID& classID()       { return CLID_LdGlastData; }

    const CsIData*    getCsIData() const    { return m_csiData; }
    const SiData*     getSiData() const{ return  m_siData; }
    const IVetoData*  getVetoData() const  { return m_vetoData; }

    const MCTruth*    getMCTruth()const { return m_truth;}
    virtual void printOn(std::ostream&)const;


    void setMCTruth(const MCTruth* t){ m_truth =t;}
private:
    const CsIData* m_csiData;
    const SiData* m_siData;
    const IVetoData* m_vetoData;
    const MCTruth*  m_truth;
};

void LdGlastData::printOn(std::ostream& out)const {
    getVetoData()->printOn(out);
    getSiData()->printOn(out);
    getCsIData()->printOn(out);
}
#endif

