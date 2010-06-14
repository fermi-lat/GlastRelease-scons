// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastHitsEvt_TrackerHits_v0_
#define _H_GlastHitsEvt_TrackerHits_v0_

// includes
#include "Gaudi/Kernel/DataObject.h"
#include "GlastEvent/Hits/SiLayerHits.h"

    
    // declare the GAUDI class-id for this class
    extern const CLID&   CLID_TrackerHits;
    
    // class TrackerHits
    /*! 
	   This class contains SiLayerHits, an array each for x- and y-orientations
	   in a single GLAST tower.
     */
    template <class Container>
        class   TrackerHits : public DataObject, 
        public HitsBase<Container> {
    public:
        //! typedef LayerList as a list of SiLayer objects
        typedef SiLayerHits<Container::SiLayerType> SiLayerType;
        typedef std::vector< SiLayerType* >         LayerList;
        
        // default constructor
        TrackerHits () : HitsBase<Container> () {}
        
        // constructor
        TrackerHits ( Container* d ) : HitsBase<Container> (d) {}
        
        // destructor
        virtual ~TrackerHits () {}  // superclasses handle all of the deletions
        
        // GAUDI identification methods
        static const CLID&    classID ()  { return CLID_TrackerHits; }
        virtual const CLID&   clID ()     { return CLID_TrackerHits; }
        
        //! Access the list of layers of x-oriented-silicon
        const LayerList&   TrackerXLayers () const { return m_xlayers; }
        
        //! Access the list of layers of y-oriented-silicon
        const LayerList&    TrackerYLayers () const { return m_ylayers; }
        
        //! Add a layer to the layers of silicon
        void    AddTrackerXLayer ( SiLayerType* t ) { m_xlayers.push_back (t); }
        
        //! Add a layer to the layers of y-silicon
        void    AddTrackerYLayer ( SiLayerType* t ) { m_ylayers.push_back (t); }
        
    private:
        
        // private data members
        
        // The list of x-oriented silicon layers in the tracker
        LayerList   m_xlayers;
        
        // The list of y-oriented silicon layers in the tracker
        LayerList   m_ylayers;
    };
    


#endif // _H_GlastHitsEvt_TrackerHits_v0_