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
#include "TdSiData.h"
#include "TdCsIData.h"
#include "TdVetoData.h"

#include "GlastEvent/data/TdGlastData.h"

#include <iostream>
#include <vector>

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
    LdGlastData ( TdCsIData* csi, TdSiData* si, TdVetoData* vd) {
        m_csiData = csi;
        m_siData = si;
        m_vetoData = vd;
    }

    virtual const CLID& clID() const   { return LdGlastData::classID(); }
    static const CLID& classID()       { return CLID_LdGlastData; }
/*
    CsIData*   getCsIData() const    { return m_csiData; }
    SiData*    getSiData() const{ return  m_siData; }
    IVetoData*  getVetoData() const  { return m_vetoData; }

 */
    const CsIData*   getCsIData() const    { return m_csiData; }
    const SiData*    getSiData() const{ return  m_siData; }
    const IVetoData*  getVetoData() const  { return m_vetoData; }
   


private:
    TdCsIData* m_csiData;
    TdSiData* m_siData;
    TdVetoData* m_vetoData;
};

#endif

