// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastEvent_EventHits
#define _H_GlastEvent_EventHits

// includes
#include "Gaudi/Kernel/DataObject.h"
#include "GlastEvent/Hits/GlastHits.h"

namespace GlastEvent { // declare the namespace for the glast raw event

    // declare the GAUDI class-id for this class
    extern const CLID&   CLID_EventHits;
    
    // class EventHits
    //! Contains a GlastHits object.
    //
    class   EventHits : public DataObject {
    public:
        /// default constructor
        EventHits () {}
        
        /// destructor
        virtual ~EventHits () {}
        
        // *** constant/type declarations ***
        
        /// GAUDI identification methods
        static const CLID&    classID ()  { return CLID_EventHits; }
        virtual const CLID&   clID ()     { return CLID_EventHits; }
        
        // the Glast containment class
        
        /// Access the glast instrument data
        const GlastHits*  InstrumentData () const { return m_glast; }
        
        /// Assign the instrument data to this event
        void    AddInstrumentData ( GlastHits* g ) { m_glast = g; }
        
    private:
        
        // private data members
        
        // The whole glast hits data underneath this object
        GlastHits*    m_glast;
        
    };

}   // namespace GlastHitsEvt


#endif // _H_GlastHitsEvt_