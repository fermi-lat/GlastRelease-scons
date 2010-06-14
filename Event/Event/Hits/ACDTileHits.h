// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastHitsEvt_ACDTileHits_v0_
#define _H_GlastHitsEvt_ACDTileHits_v0_

// includes
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/Hits/HitsBase.h"


    // declare the GAUDI class-id for this class
    extern const CLID&   CLID_ACDTileHits;
    
    // class ACDTileHits
    //!  This class represents the data associated with a single silicon strip hit in the tracker. 
    /*! 
	 Note that this is a template class, which
     will must be instanciated with an ACDTile representation. (Currently instrument/Scintillator.)
    
     The template argument class must support the following methods:
    
      - ACCESS
         - float   energy () const  [ Access the energy in the tile ]
             
      - SET
         - void    energy ( float ) [ Set the energy in the tile ]
    */
    template <class Container>
    class   ACDTileHits :   public ContainedObject,         // for Gaudi
                            public HitsBase <Container> {   // for our containment management
    public:
        // default constructor (will allocate automatically)
        ACDTileHits () : HitsBase<Container> () {}
        
        // constructor
        ACDTileHits (Container* s) : HitsBase<Container> (s) {}
        
        // destructor
        virtual ~ACDTileHits () {}
        
        // GAUDI identification methods
        static const CLID&    classID ()  { return CLID_ACDTileHits; }
        virtual const CLID&   clID ()     { return CLID_ACDTileHits; }

        //! Access the energy deposited in the tile (in GeV)
        //! A negative value indicates an error.
        float  Energy () const { return haveObj() ? getContainer()->energy() : -1.; }

        //! set the various parameters associated with this class
        void   SetEnergy ( float e ) { if (haveObj()) getContainer()->energy(e); }

    private:
        
        // private data members
        
    };



#endif // _H_GlastHitsEvt_ACDTileHitsData_v0_