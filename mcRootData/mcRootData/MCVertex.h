// McVertex.h
//
// The McVertex class is intended to contain monte carlo particle position data for the 
// GLAST mc.  Designed initially to meet the simple needs of the GLAST balloon test.
//
// May 3 2001 Daniel Flath - Creation

#ifndef MCVERTEX_H
#define MCVERTEX_H

#include "TObject.h"
#include "TLorentzVector.h"

/*! GLAST Monte Carlo vertex class.  Designed originally for 2001 balloon flight.
 *  Gives position information for MCParticle which keeps a pointer to a vertex.
 *  
 *  \author Daniel Flath
 *  
 *  \b Revisions:
 *  - 03 May 2001   Daniel Flath    Creation, "get" functions added
 *  - 11 May 2001   Daniel Flath    "set" functions added
 *  - 14 May 2001   Daniel Flath    DOxygen style comments added
 */
class MCVertex: public TObject {

  private:
    //! position/time 4-vector
    TLorentzVector *m_pPos;

  public:
    ////////////////////// construction/destruction: ///////////////////////

    //! default constructor. DO NOT USE -- see details.
    MCVertex();

    //! constructor
    MCVertex(TLorentzVector *pos);

    //! constructor
    MCVertex(Double_t posX, Double_t posY, Double_t posZ, Double_t time);

    //! destructor
    ~MCVertex();

    //////////////////////   data access functions:  ///////////////////////

    //! get position x-component
    Double_t getPosX();

    //! get position y-component
    Double_t getPosY();

    //! get position z-component
    Double_t getPosZ();

    //! get time
    Double_t getTime();

    //! get position/time 4-vector
    TLorentzVector *getPosLV();


    //! set position x-component
    void setPosX(Double_t posX);

    //! set position y-component
    void setPosY(Double_t posY);

    //! set position z-component
    void setPosZ(Double_t posZ);

    //! set time
    void setTime(Double_t time);

    ClassDef(MCVertex,1)
};

#endif