// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastEvent_CalorimeterHits_v0_
#define _H_GlastEvent_CalorimeterHits_v0_

// includes
#include "Gaudi/Kernel/DataObject.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Hits/CalorimeterLogHits.h"

    
    // declare the GAUDI class-id for this class
    extern const CLID&   CLID_CalorimeterHits;
    
    // class CalorimeterHits
    //!  Wraps a class containing a list of calorimeter log hit objects. 
	/*!  
	    The wrapped class must typedef "LogType". (Currently instrument/Calorimeter).

    */
    template <class Container>
        class   CalorimeterHits : public ::DataObject,            // for Gaudi
        public HitsBase<Container> {  // for our own management
    public:
        //! typedef the list of calorimeter layers
        typedef ObjectVector< CalorimeterLogHits<Container::LogType> >  LogList;
        
        // default constructor
        CalorimeterHits () : HitsBase<Container> () {}
        
        // constructor
        CalorimeterHits ( Container* d );
        
        // destructor
        virtual ~CalorimeterHits () {} // ObjectVector auto-destroys!
        
        // GAUDI identification methods
        static const CLID&    classID ()  { return CLID_CalorimeterHits; }
        virtual const CLID&   clID ()     { return CLID_CalorimeterHits; }
        
        //! Access the list of calorimeter layers
        const LogList&    CalLogs () const { return m_layers; }
        
    private:
        
        // private data members
        
        // the list of calorimeter layers
        LogList   m_logs;
    };

    // inline methods
    
    // constructor - assigns to the object vector 
    template <class Container>
        CalorimeterHits<Container>::CalorimeterHits ( Container* d )
        : HitsBase <Container> (d)
    {
        // now run through and create a list of CalorimeterLogHits containers
        // for the actual log data we contain.
        Container::const_iterator it = d->begin(); 
        for (; it != d->end(); ++it) {
            
            // access the object
            m_logs->push_back (*it);
        }
    }



#endif // _H_GlastHitsEvt_CalorimeterHits_v0_