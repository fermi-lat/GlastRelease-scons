// RMcParticle.h
//
// The RMcParticle class is intended to contain monte carlo particle data for the 
// GLAST MC.  Designed initially to meet the simple needs of the GLAST balloon flight sim.
//
// May 3 2001 Daniel Flath - Creation

#ifndef RMCPARTICLE_H
#define RMCPARTICLE_H

#include "TObject.h"
#include "TObjArray.h"
#include "TLorentzVector.h"

#include "RMcVertex.h"


/*! GLAST Monte Carlo particle class.  Designed originally for 2001 balloon flight.
 *  Inheritance information given for parentage only.
 *
 *  \author Daniel Flath
 *  
 *  \b Revisions:
 *  - 03 May 2001   Daniel Flath    Creation, "get" functions added
 *  - 11 May 2001   Daniel Flath    "set" functions added
 *  - 14 May 2001   Daniel Flath    DOxygen style comments added
 *  - 06 Jun 2001   Daniel Flath    Final rewrite for integration with ROOTWriter
 *  - 14 Jun 2001   Daniel Flath    DOxygen style comments updated for checkin with ROOTWriter
 */
class RMcParticle: public TObject {
private:
    //! particle id number
    Int_t m_nPDGId;

    //! initial position
    RMcVertex *m_pVert;

    //! parent particle
	RMcParticle *m_pParent;

    //! list of child particles
    /*! Ensure that if list is modified, that m_nChildCount is updated accordingly
     */
    TObjArray *m_pChildList;

    //! number of child particles in list
    /*! Must be updated to reflect modifications to the m_pChildList
     */
    Int_t m_nChildCount;

public:
    ////////////////////// construction/destruction: ///////////////////////

    //! default constructor. DO NOT USE -- see details.
    RMcParticle();

    //! constructor
    RMcParticle(Int_t id, RMcParticle *pParent, RMcVertex *pVertex);

    //! destructor
    ~RMcParticle();

    //////////////////////   data access functions:  ///////////////////////

    //! get pointer to associated vertex class instance
    inline RMcVertex *getVertex()    const { return m_pVert; };

    //! get parent particle
    inline RMcParticle *getParent()  const { return m_pParent; };

    //! get particle id number
    inline Int_t getPDGId()         const { return m_nPDGId; };

    //! get number of children
    inline Int_t getChildCount()    const { return m_nChildCount; };

    //! get child from list at index
    RMcParticle *getChild(Int_t nIndex) const;    

    //! add a child to list
    void addChild(RMcParticle *pChild);


    ClassDef(RMcParticle,1)
};

#endif