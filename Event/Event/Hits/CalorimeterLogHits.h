// $Id$
//
//  Original Author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastEvent_CalorimeterLogHits_v0_
#define _H_GlastEvent_CalorimeterLogHits_v0_

// includes
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/Hits/HitsBase.h"


    // declare the GAUDI class-id for this class
    extern const CLID&   CLID_CalorimeterLogHits;
    
    // class CalorimeterLogHits
    //!  Provide access to the Log hits for a single calorimeter log.
	/*! The wrapped class Container currently must provide:
	     - int layer() const
		 - int layer(int i)
		 - float Lresp() const;
		 - float Rresp() const;
		 - float energy() const;
		 - void loadXtal( float E, float plus, float minus );
     */
    template <class Container>
        class   CalorimeterLogHits : public ContainedObject,
        public HitsBase<Container> {
    public:
        // default constructor
        CalorimeterLogHits () : HitsBase<Container> () {}
        
        // constructor
        CalorimeterLogHits (Container* l) : HitsBase<Container> (c) {}
        
        // destructor
        virtual ~CalorimeterLogHits () {}
        
        // *** constants/enumerations ***
        
        // enumerator to distinguish between the two ends of the log
        typedef enum {
            Plus = 0,
                Minus
        } LogEnd;
        
        // GAUDI identification methods
        static const CLID&    classID ()  { return CLID_CalorimeterLogHits; }
        virtual const CLID&   clID ()     { return CLID_CalorimeterLogHits; }
        
        // *** Access methods ***
        
        //! Access the layer within the Calorimeter for a the log
        //!  negative value indicates error.
        int     Layer () const { return hasObj() ? getContainer()->layer() : -1; }
        
        //! Access the index (left->right position) of the log
        //!  negative value indicates error.
        int     Index () const { return hasObj() ? getContainer()->index() : -1; }
        
        //! Access the plus or minus response for the log
        float   EndResp ( LogEnd e ) const 
        { return hasObj() ? ((e == Plus) ? getContainer->Lresp() : getContainer()->Rresp()) : -1; }
        
        //! Access the total energy deposition (GeV)
        float   TotE () const { return hasObj() ? getContainer()->energy() : -1.; }
        
        // *** Set methods ***
        
        //! set the layer for this log
        void    Layer ( int i ) { if (hasObj()) getContainer()->layer (i); }
        
        //! set the index for this log
        void    Index ( int i ) { if (hasObj()) getContainer()->index (i); }
        
        //! set response for this log - total deposition, plus, minus responses
        void    SetResponse ( float E, float plus, float minus ) 
        { if (hasObj()) getContainer()->loadXtal (E, plus, minus); }
        
    private:
        
        // private data members
        
    };
    

#endif // _H_GlastHitsEvt_CalorimeterLogHits_v0_