// TdVetoData.h: interface for the TdVetoData class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TDVETO_DATA_H
#define TDVETO_DATA_H 1

#include "data/IVetoData.h"

#include "instrument/Glast.h"
#include "geometry/Point.h"

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/CellID.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"


// forward declarations
class Scintillator;

extern const CLID& CLID_TdVetoData;

//Original Author: Ian Gable (03/02/2001)

/*! Implements the IVetoData interface. Introduced primarily for the purpose
    of getting a functional aoRecon
*/
class TdVetoData : public IVetoData, public DataObject
{
public:

    virtual const CLID& clID() const   { return TdVetoData::classID(); }
    static const CLID& classID()       { return CLID_TdVetoData; }

    TdVetoData ();
    ~TdVetoData ();
    
    TdVetoData::const_iterator begin () const
    {
        return tileList.begin();
    }
    
    TdVetoData::const_iterator end () const
    {
        return tileList.end();
    }

    //! set from GlastDetector objects
    void load (const Scintillator& tile);
    
    //! formatted dump
    void printOn (std::ostream& cout) const;
    
    
    void clear ();

    int count()const;

private:
    
    //        storage is a list of Tiles, defined in interface
    VetoList tileList;
    
private:
};
 
#endif

