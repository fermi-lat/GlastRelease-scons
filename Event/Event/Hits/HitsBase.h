// $Id$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastEvent_HitsBase_
#define _H_GlastEvent_HitsBase_

// declare everything in the GlastEvent namespace

    // class HitsBase
    /*!  This is a wrapper around a class (Container) that has hit data.

          The Container must have a public default constructor and destructor.
     */
    template <class Container> 
        class   HitsBase {
    public:
        
        // *** ctor/dtor ***
        
        /// default constructor (allocates a container on true)
        HitsBase ( bool alloc = true )
            : m_owns (true), m_container ((alloc) ? new Container : 0) {}
        
        /// constructore - assigns a container
        HitsBase ( Container* c, bool owns = true ) 
            : m_owns (owns), m_container (c) {}
        
        /// destructor - deletes the container according to the owns
        virtual ~HitsBase () { if (m_owns) delete m_container; m_container = 0; }
        
        // *** access methods ***
        
        /// do we own the container? 
        bool    ownsContainer () const { return m_owns; }
        
        /// do we have a valid object? (non null)
        bool    haveObj () const { return m_container != 0; }
        
        /// access the container (const)
        const Container*    operator () (void) const { return getContainer(); }
        
        // access the container (non-const)
        Container*  operator () (void) { return getContainer(); }
        
    protected:
        
        // *** protected access ***
        
        /// access the container (internal const)
        const Container*    getContainer () const { return m_container; }
        
        /// access the container (internal non-const)
        Container*  getContainer () { return m_container; }
        
    private:
        
        // *** data members ***
        
        // do we own the container? 
        bool    m_owns;
        
        // the container object itself
        Container*  m_container;
    };
    
    // inline functions
#endif  // _H_GlastEvent_HitsBase_