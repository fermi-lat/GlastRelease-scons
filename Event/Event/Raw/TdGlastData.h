// TdGlastData.h: interface for the TdGlastData class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TdGlastData_H
#define TdGlastData_H 1

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

#include "instrument/CsIDetector.h"
#include "data/GlastData.h"


// forward declarations
#include "GlastEvent/Raw/TdSiData.h"
#include "GlastEvent/Raw/TdCsIData.h"

#include <iostream>
#include <vector>

extern const CLID& CLID_TdGlastData;

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


class TdGlastData : virtual public GlastData , virtual public DataObject {

public:
    TdGlastData ( TdCsIData* csi, TdSiData* si) {
        m_csiData = csi;
        m_siData =si;
    }

    virtual const CLID& clID() const   { return TdGlastData::classID(); }
    static const CLID& classID()       { return CLID_TdGlastData; }


    const CsIData*   getCsIData() const    { return m_csiData; }
    const SiData*    getSiData() const{ return  m_siData; }
    const IVetoData*  getVetoData() const  { return NULL; }


private:
    TdCsIData* m_csiData;
    TdSiData* m_siData;
    // Note that TdVetoData will need to be added here at 
    // some point.
};

#endif
