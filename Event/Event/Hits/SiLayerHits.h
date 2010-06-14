// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastEvent_SiLayerHits_v0_
#define _H_GlastEvent_SiLayerHits_v0_

// includes
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Hits/SiStripHits.h"
#include "GlastEvent/Hits/HitsBase.h"

    
    // declare the GAUDI class-id for this class
    extern const CLID&   CLID_SiLayerHits;
    
    
    // class SiLayerHits
    //!  Wraps a class representing a single Si hit.
    /*!
	     The class must implement the methods
		 - int getId() const
		 - bool isXplane() const
		 - int size() const

        And must typedef SiStripType. Currently instrument/SiDetector satisfies this, 
		and SiStripType is then SiDetector::Strip.
     */
    template <class Container>
        class   SiLayerHits : public DataObject,
        public HitsBase<Container> {
    public:
        //! HitList - typedef the list of hits for a single silicon layer
        typedef ObjectVector < SiStripHits<Container::SiStripType> >    StripList;
        
        // default constructor
        SiLayerHits () : HitsBase<Container> () {}
        
        // constructor - assign object
        SiLayerHits (Container* d) : HitsBase<Container> (d) {}
        
        // destructor
        virtual ~SiLayerHits () {}  // super classes handle deletion
        
        
        // *** type declarations/indexes ***
        
        //! declare a type to handle the indexing of the layers
        typedef struct {
            int     index;  // index of the layer from the bottom
            bool    isX;    // is this an X layer ?
        } IdType;
        
        // GAUDI identification methods
        static const CLID&    classID ()  { return CLID_SiLayerHits; }
        virtual const CLID&   clID ()     { return CLID_SiLayerHits; }
        
        // *** access methods ***
        
        //! Access the vertical index for this layer (bottom to top)
        int     Index () const { return hasObj() ? getContainer()->getId() : Container::invalid_id(); }
        
        //! Is this an X-oriented plane (strips || to x-axis, measuring Y scalar) 
        bool    IsXPlane () const { return hasObj() ? getContainer()->isXplane() : 0; }
        
        //! Access the number of hits seen by one controller
        int   NHits () const { return hasObj() ? getContainer()->size() : 0; }
        
        //! Access the hits associated with this layer
        const StripList&  Hits () const { return m_hits; }
        
    private:
        
        // private data members
        
        // The hit list which wraps the hits within the SiLayerHits object
        StripList     m_hits;
    };
    


#endif // _H_GlastHitsEvt_SiLayerHits_v0_