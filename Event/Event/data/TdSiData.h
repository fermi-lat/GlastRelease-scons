// $Id$

#ifndef RECOSI_DATA_H
#define RECOSI_DATA_H 1

//#include "instrument/Glast.h"
#include "geometry/Point.h"
#include "idents/ModuleId.h"
#include "data/SiData.h"

#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/CellID.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"


#include <iostream>
#include <vector>

extern const CLID& CLID_TdSiData;

class TdSiData : virtual public SiData , virtual public DataObject {
        public:
    
    class Strip {
        private:
            friend class TdSiData;
        
        public:
            Strip () {}
        
        protected:
            Strip (Point p, Id m, unsigned int n, unsigned int t = 0)
                : stripIndex(n), stripType(t), pos(p), module(m)
            {}
        
            unsigned int stripIndex;
            unsigned int stripType; 
            Point pos;
            Id module;
    };

    virtual const CLID& clID() const   { return TdSiData::classID(); }
    static const CLID& classID()       { return CLID_TdSiData; }

        
    TdSiData ();
    
    TdSiData (unsigned int n);
        
    virtual ~TdSiData ();
    
    //	access to module id if desired
    virtual const SiData::Id& moduleId (enum SiData::Axis a, 
        unsigned int tray, 
        unsigned int n) const;
    
    
    
    //-------------------- access functions for data analysis------------------------
    void clear ();
    
    //	number of hit strips in the given tray (tray=0..16)
    int nHits (enum SiData::Axis a, int tray) const;
    
    //	center of the given strip, where  n=0..nHits(a,tray,n)-1
    Point hit (enum SiData::Axis a, unsigned int tray, unsigned int n) const;
    
    //	strip Id
    unsigned int hitId (enum SiData::Axis a, 
        unsigned int tray, unsigned int n) const;
    
    //	Strip Type (1 = noise, 0 = real)
    //unsigned int hitType (enum SiData::Axis a, 
    //    unsigned int tray, unsigned int n) const;
    
    //	total number of Si strips hit
    //	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int totalHits () const;
    
    //	largest number of hits on an individual readout
    //	controller
    //	(determines readout time)
    //	----- I/O Functions -----------------------
    //	total number of Si strips hit
    //int maxControllerHits () const;
    
    void readData (std::istream& in);
    
    //	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //	Data I/O
    void writeData (std::ostream& out);
    
    //	---------------------------------------------------------
    void printOn (std::ostream& cout = std::cout) const;
    
protected:
    
    int m_total_hits;
    int m_controller_max;
//    class Strip;
    //	storage is a vector (by layer number) of lists of hit strips
    std::vector< std::vector< class Strip >* > xhitList;
    
    //	storage is a vector (by layer number) of lists of hit strips
    std::vector< std::vector< class Strip >* > yhitList;
    
};

#endif
