// RMcVertex.h
//
// The RMcVertex class is intended to contain monte carlo particle position data for the 
// GLAST MC.  Designed initially to meet the simple needs of the GLAST balloon flight sim.
//
// May 3 2001 Daniel Flath - Creation

#ifndef RMCVERTEX_H
#define RMCVERTEX_H

#include "TObject.h"
#include "TLorentzVector.h"

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
 */
class RMcVertex: public TObject {

  private:
    //! momentum/energy 4-vector
    TLorentzVector *m_pMomI, *m_pMomF;

    //! position/time 4-vector
    TLorentzVector *m_pPosI, *m_pPosF;

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


    ClassDef(RMcVertex,1)
};

#endif