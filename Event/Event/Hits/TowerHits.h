// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastEvent_TowerHits_v0_
#define _H_GlastEvent_TowerHits_v0_

// includes
#include "Gaudi/Kernel/DataObject.h"
#include "GlastEvent/Hits/ACDTileHits.h"
#include "GlastEvent/Hits/CalorimeterHits.h"
#include "GlastEvent/Hits/TrackerHits.h"
#include "GlastEvent/Hits/HitsBase.h"

    
    // declare the GAUDI class-id for this class
    extern const CLID&   CLID_TowerHits;
    
    // class TowerHits
    /*!  This contains MC hits for a single tower. It has points to associated
	     TrackerHits and CalorimeterHits objects.
    */
    template <class Container>
        class   TowerHits : public DataObject,
        public HitsBase<Container> {
    public:
        // default constructor
        TowerHits () : HitsBase <Container> () {}
        
        // constructor
        TowerHits ( Container* d ) : HitsBase <Container> (d) {}
        
        // destructor
        virtual ~TowerHits () {}  // super classes take care of deletion
        
        // *** type/constant declarations ***
        
        //! declare a type to hold the tracker information
        typedef TrackerHits<Container::TrackerType> TrackerHitsType;
        
        //! declare a type to hold the calorimeter information
        typedef CalorimeterHits<Container::CalorimeterType> CalorimeterHitsType;
        
        // GAUDI identification methods
        static const CLID&    classID ()  { return CLID_TowerHits; }
        virtual const CLID&   clID ()     { return CLID_TowerHits; }
        
        //! Access to the tracker data for this tower
        const TrackerHitsType*        GetTrackerHits () const { return m_tracker; }
        
        //! Access to the calorimeter data for this tower
        const CalorimeterHitsType*    GetCalorimeterHits () const { return m_calorimeter; }
        
        //! Set the calorimeter hits 
        void    SetCalorimeterHits ( CalorimeterHitsType* cal ) { m_calorimeter = cal; }
        
        //! Set the tracker hits
        void    SetTrackerHits ( TrackerHitsType* tkr ) { m_tracker = tkr; }
        
        //! Index of this tower
        int     Index () const { return hasObj() ? getContainer()->getId() : -1; }
        
    private:
        
        // private data members
        
        // The raw data for the tracker
        TrackerHits*         m_tracker;
        
        // The raw data for the calorimeter
        CalorimeterHits*     m_calorimeter;
        
    };
    


#endif // _H_GlastHitsEvt_TowerHits_v0_