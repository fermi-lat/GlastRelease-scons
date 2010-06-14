// McParticle.h
//
// The McPart class is intended to contain monte carlo particle data for the 
// GLAST mc.  Designed initially to meet the simple needs of the GLAST balloon test.
//
// May 3 2001 Daniel Flath - Creation

#ifndef MCPARTICLE_H
#define MCPARTICLE_H

class MCVertex;

#include "TObject.h"
#include "TLorentzVector.h"

/*! GLAST Monte Carlo particle class.  Designed originally for 2001 balloon flight.
 *  Inheritance information given for parentage only.
 *  
 *  \author Daniel Flath
 *  
 *  \b Revisions:
 *  - 03 May 2001   Daniel Flath    Creation, "get" functions added
 *  - 11 May 2001   Daniel Flath    "set" functions added
 *  - 14 May 2001   Daniel Flath    DOxygen style comments added
 */
class MCParticle: public TObject {

  private:
    //! momentum/energy 4-vector
    TLorentzVector *m_pMom;

    //! particle id number
    Int_t m_nPDGId;

    //! initial position
    MCVertex *m_pVert;

    //! parent particle
	  MCParticle *m_pParent;

  public:
    ////////////////////// construction/destruction: ///////////////////////

    //! default constructor. DO NOT USE -- see details.
    MCParticle();

    //! constructor
    MCParticle(TLorentzVector *pMom, Int_t id, MCParticle *pParent, MCVertex *pVertex);

    //! constructor
    MCParticle(Double_t momX, Double_t momY, Double_t momZ, Double_t energy,
               Int_t id, MCParticle *pParent, MCVertex *pVertex);

    //! destructor
    ~MCParticle();

    //////////////////////   data access functions:  ///////////////////////

    //! get momentum x-component
    Double_t getMomX();

    //! get momentum y-component
    Double_t getMomY();

    //! get momentum z-component
    Double_t getMomZ();

    //! get magnitude of total momentum
    Double_t getMom();

    //! get energy
    Double_t getEnergy();

    //! get momentum/energy 4-vector
    TLorentzVector *getMomLV();

    //! get pointer to vertex class instance
    MCVertex *getVertex();

    //! get parent particle
	  MCParticle *getParent();

    //! get particle id number
    Int_t getPDGId();

    
    //! set momentum x-component
    void setMomX(Double_t momX);

    //! set momentum y-component
    void setMomY(Double_t momY);

    //! set momentum z-component
    void setMomZ(Double_t momZ);

    //! set energy
    void setEnergy(Double_t energy);

    ClassDef(MCParticle,1)
};

#endif