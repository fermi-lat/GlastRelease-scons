// RMcVertex.h
//
// The RMcVertex class is intended to contain monte carlo pa
//! position/time 4-vectorrticle position data for the 
// GLAST MC.  Designed initially to meet the simple needs of the GLAST balloon flight sim.
//
// May 3 2001 Daniel Flath - Creation

#ifndef RMCVERTEX_H
#define RMCVERTEX_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TObjArray.h"

class RMcParticle;


/*! GLAST Monte Carlo vertex class.  Designed originally for 2001 balloon flight.
*  Gives position information for McParticle which keeps a pointer to a vertex.
*  
*  \author Daniel Flath
*  
*  \b Revisions:
*  - 03 May 2001   Daniel Flath    Creation, "get" functions added
*  - 11 May 2001   Daniel Flath    "set" functions added
*  - 14 May 2001   Daniel Flath    DOxygen style comments added
*  - 06 Jun 2001   Daniel Flath    Final rewrite for integration with ROOTWriter
*  - 14 Jun 2001   Daniel Flath    DOxygen style comments updated for checkin with ROOTWriter
*  - 14 Aug 2001   Ian Gable       Rehack to match the TDS versions of these classes
*/
class RMcVertex: public TObject {
	
private:
    //! momentum/energy 4-vector
    TLorentzVector *m_pMomI, *m_pMomF;
	
    TLorentzVector *m_pPosI, *m_pPosF;
	
	Double_t m_timeOfFlight;
	//!Array of Daughter particles
	TObjArray* m_mcDaughterParticles;
	//! Mother particles
	RMcParticle* m_motherParticle;
    //!Particle associated with this Vertex
    RMcParticle* m_mainParticle;
	
public:
    ////////////////////// construction/destruction: ///////////////////////
	
    //! default constructor. DO NOT USE -- see details.
    RMcVertex();
	
    //! constructor
    RMcVertex(TLorentzVector *pPosI, TLorentzVector *pPosF, 
		TLorentzVector *pMomI, TLorentzVector *pMomF);
	
    //! destructor
    ~RMcVertex();
	
    //////////////////////   data access functions:  ///////////////////////
    
    //! get initial position/time 4-vector
    const TLorentzVector *getPosInitial() const { return m_pPosI; };
	
    //! get final position/time 4-vector
    const TLorentzVector *getPosFinal()   const { return m_pPosF; };
	
    //! get initial momentum/energy 4-vector
    const TLorentzVector *getMomInitial() const { return m_pMomI; };
	
    //! get final momentum/energy 4-vector
    const TLorentzVector *getMomFinal()   const { return m_pMomF; };
	
    //! set the Initial Position
	void setPosInitial(TLorentzVector* pPosI) { m_pPosI = pPosI; };

	//! set the Final Position
	void setPosFinal(TLorentzVector* pPosF) { m_pPosF = pPosF; };

	//! set the Initial Momentum
	void setMomInitial(TLorentzVector* pMomI) {m_pMomI = pMomI; };

	//! set the Final Momentum
	void setMomFinal(TLorentzVector* pMomF) { m_pMomF = pMomF; };
	
	
	/// retrieve time of flight
    Double_t timeOfFlight () const { return m_timeOfFlight;};

    /// update time of flight
    void setTimeOfFlight (double value) {m_timeOfFlight = value;};
	
    /// Retrieve pointer to the pair particle 
	RMcParticle* mcParticle();
    /// Update pointer to the pair particle (by a C++ pointer)
    void setMcParticle( RMcParticle* value );
	
    /// Retrieve pointer to mother particle (const )
    RMcParticle* motherMcParticle();
	
    /// Update pointer to mother particle (by a C++ pointer or a smart reference)
    void setMotherMcParticle( RMcParticle* value );
	
    /// Retrieve pointer to vector of daughter particles (const or non-const)
    const TObjArray* daughterMcParticles() const;
	
    /// add Paticles to the TObjArray
    void addDaughterMcParticle( RMcParticle* particle );
	
	
	
    ClassDef(RMcVertex,1)
};

#endif