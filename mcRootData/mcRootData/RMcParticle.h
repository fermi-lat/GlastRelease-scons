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

class RMcVertex;

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
 *  - 23 Aug 2001   Ian Gable       Rehacked to work with the TDS classes
 */
class RMcParticle: public TObject {

private:
    //! particle id number
    Int_t m_particleId;

    //! initial position
    RMcVertex *m_mcVertex;

    //! Pramary particle indicator
    Bool_t m_primary;


public:
    ////////////////////// construction/destruction: ///////////////////////

    //! default constructor. DO NOT USE -- see details.
    RMcParticle();

    //! constructor
    RMcParticle(Int_t id, RMcVertex *pVertex);

    //! destructor
    ~RMcParticle();

    //////////////////////   data access functions:  ///////////////////////

    //! get particle id number
    Int_t particleId();

    /// Update particle identification
    void setParticleId( Int_t value );

    /// Retrieve whether this is a primary particle
    Bool_t primaryParticle();

    // Set whether this is a primary particle
    void setPrimaryParticleFlag( Bool_t value );


	 /// Retrieve pointer to the vertex
    RMcVertex* mcVertex();

    /// Update pointer to origin vertex (by a C++ pointer or a smart reference)
    void setMcVertex( RMcVertex* value );

    ClassDef(RMcParticle,1)
};

#endif