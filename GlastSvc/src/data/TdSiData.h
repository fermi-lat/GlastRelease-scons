#ifndef RECOSI_DATA_H
#define RECOSI_DATA_H 1

#include "instrument/Glast.h"
#include "geometry/Point.h"
#include "idents/ModuleId.h"
#include "data/SiData.h"

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/CellID.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"


#include <iostream>
#include <vector>

extern const CLID& CLID_TdSiData;

class TdSiData : virtual public SiData , virtual public DataObject {
public:

    virtual const CLID& clID() const   { return TdSiData::classID(); }
    static const CLID& classID()       { return CLID_TdSiData; }

        
    TdSiData ();
    
    TdSiData (unsigned int n);
    
    //copy constructor
    TdSiData (TdSiData*,unsigned int n);
    
    virtual ~TdSiData ();
    
    //  access to module id if desired
    virtual const SiData::Id& moduleId (enum SiData::Axis a, 
        unsigned int tray, 
        unsigned int n) const;
    
    
    void load (const SiDetector& strips, idents::ModuleId module);
    
    //-------------------- access functions for data analysis------------------------
    void clear ();

    //  number of trays
    int nTrays(enum SiData::Axis a) const;
    
    //  number of hit strips in the given tray (tray=0..16)
    int nHits (enum SiData::Axis a, int tray) const;
    
    //  center of the given strip, where  n=0..nHits(a,tray,n)-1
    Point hit (enum SiData::Axis a, unsigned int tray, unsigned int n) const;
    
    //  strip Id
    unsigned int hitId (enum SiData::Axis a, 
        unsigned int tray, unsigned int n) const;
    
    //  Strip Type (1 = noise, 0 = real)
    //unsigned int hitType (enum SiData::Axis a, 
    //    unsigned int tray, unsigned int n) const;
    
    //  total number of Si strips hit
    int totalHits () const;
    
    
    //  ---------------------------------------------------------
    void printOn (std::ostream& cout = std::cout) const;
    
private:
    
    int m_total_hits;
    int m_controller_max;
    class Strip;
    //  storage is a vector (by layer number) of lists of hit strips
    std::vector< std::vector< class Strip >* > xhitList;
    
    //  storage is a vector (by layer number) of lists of hit strips
    std::vector< std::vector< class Strip >* > yhitList;
    
};

#endif
