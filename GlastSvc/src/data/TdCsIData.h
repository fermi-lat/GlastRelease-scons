//old school version
#ifndef TDCSI_DATA_H
#define TDCSI_DATA_H 1


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
#include "GlastEvent/Utilities/CellID.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"

#include "instrument/CsIDetector.h"

#include "data/CsIData.h"

class CsIDetector;

extern const CLID& CLID_TdCsIData;

//! Representation of the CsI Calorimeter Data
/*! \class TdCsIData
*/

class TdCsIData : virtual public DataObject , virtual public CsIData
{

public:
    virtual const CLID& clID() const   { return TdCsIData::classID(); }
    static const CLID& classID()       { return CLID_TdCsIData; }
    
    
public:
/*! Inner class used to store data on the individual crystals
    */
    class Xtal 
    {
        
        
    private:
        friend class TdCsIData;
    public:
        
        Xtal (Point p, float e, idents::ModuleId m, int i, float Left, float Right, std::vector<double> diodes_energy)
            : pos(p), energy(e), id(i), Lresp(Left), Rresp(Right), module(m), Diodes_Energy(diodes_energy)
        {}
        Xtal (Point p, float e, idents::ModuleId m, int i, float Left, float Right)
            : pos(p), energy(e), id(i), Lresp(Left), Rresp(Right), module(m), Diodes_Energy(4,0.)
        {}
        
        
        
    private:
        
        //! Data Members for Class Attributes
        Point pos;
        float energy;
        int id;
        float Lresp;
        float Rresp;
        idents::ModuleId module;
        std::vector<double> Diodes_Energy;
    };
    
    
    
    public:
        TdCsIData(){};
        TdCsIData (int numLayers);
        void copyUp (TdCsIData* copy, int numLayers);
        
         ~TdCsIData ();
        
        //!     number of hit towers in the layer layer=0..
        int nHits (unsigned int layer) const;
        
        //!     energy recorded in the ith tower
        float energy (unsigned int layer, unsigned int n) const;
        
        //!     position of the center of the crystal
        Point xtalPos (unsigned int layer, unsigned int n) const;
        
        //!     access to module id if desired
        idents::ModuleId moduleId (unsigned int layer, unsigned int n) const;
        
        
        //!     access to xtal  id if desired
        idents::XtalId xtalId (unsigned int layer, unsigned int n) const;
        
        float Lresp (unsigned int layer, unsigned int n) const;
        
        //!     Left and Right response of this xtal
        float Rresp (unsigned int layer, unsigned int n) const;
        
        //!   energy deposited in each of 4 diodes (left-big, right-big, left-small,right-small)
        const std::vector<double>& Diodes_Energy(unsigned int layer, unsigned int n) const; 
        
        //!load local array of xtal objects from a CsIDetector
        void load (const CsIDetector& xtal, idents::ModuleId tower);
        
        //     const std::vector<double>& Diodes_Energy(unsigned int layer, unsigned int n);
        
                
        void clear ();
        
        
             
        void printOn (std::ostream& cout) const ;
        
    private:
        
        
        //!     data: id followed by left, right response (keV units)
        std::map< idents::XtalId, std::pair< int, int > > m_xtals;
        
        //!     storage is a list of lists (deprecated)
        std::vector< std::vector< class Xtal >* > calorList;
        
       
        
};

        
#endif

