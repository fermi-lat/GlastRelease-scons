// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastHitsEvt_GlastHits_v0_
#define _H_GlastHitsEvt_GlastHits_v0_

// includes
#include "Gaudi/Kernel/DataObject.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Hits/ACDhit.h"

    
    // declare the GAUDI class-id for this class
    const static CLID   CLID_GlastHits = 201;
    
    // class GlastHits
    //! This is a DataObject that contains TowerHits and ACDTileHits.
    /*! 
	   It is instantiated with a Container that knows about TowerType and ADCTileType? 
	   Not clear why it would inherit from HitsBase, since it does not contain hits itself.

     */
    class   GlastHits : public DataObject {
    public:
        
        // default constructor
        GlastHits ()  {}
                
        // destructor
        virtual ~GlastHits () {}    // superclasses handle all of the other methods
        
        // GAUDI identification methods
        static const CLID&    classID ()  { return CLID_GlastHits; }
        virtual const CLID&   clID ()     { return CLID_GlastHits; }
                
        //! Access the list of ACDTile objects
        const ACDhitVector&  ACDTiles () const { return m_tiles; }
                
        //! Add an ACDTileHits to the list of ACDTiles
        void    AddACDTile ( ACDhit* t ) { m_tiles.push_back(t); }
        
        
    private:
        
        // private data members
                
        // List of ACDTile objects underneath this object
        ACDhitVector m_tiles;
    };
    


#endif // _H_GlastHitsEvt_TowerData_v0_