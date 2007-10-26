#ifndef McTkrStrip_H
#define McTkrStrip_H

//include files

#include "TObject.h"
#include "TVector.h"
#include "TRef.h"

#include "commonRootData/idents/VolumeIdentifier.h"
#include "McPositionHit.h"

/*!
* \class McTkrStrip
* \author T. Burnett
* \brief Represent energy deposited (or charge) per Si Strip
*
*  Intermediate class for TkrDigi process.



<table>
<tr><th> Method <th> Return Type <th> Description <th>

<tr><td> getVolumeId
<td>const VolumeIdentifier<td> Returns the volume identifier
<tr><td> getTotalEnergy
<td>Double_t<td> Returns the total energy deposited in a volume
<tr><td> getMoment1
<td>const TVector3<td> Returns the first moment
<tr><td> getMoment2
<td>const TVector3<td> Returns the second moment
<tr><td> getMcParticleEnergy( Particle p )
<td>Double_t<td> Returns energy associated with PRIMARY, ELECTRON or POSITRON

</table>
*/

class McTkrStrip : public TObject
{   
public:
    //! Constructors
    //! Null constructor
    McTkrStrip();
    
    //! constructor with plane id.
    McTkrStrip(const VolumeIdentifier& id, UInt_t strip, Double_t e=0,
         Bool_t noise=0, Double_t deltaX=0, Double_t deltaY=0, TRefArray* hitList = 0);
    //! Destructor
    virtual ~McTkrStrip();
    
    //! Initialize
    void initialize(const VolumeIdentifier& id, UInt_t strip, Double_t e,
         Bool_t noise, Double_t deltaX, Double_t deltaY, TRefArray* hitList);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    //! access to the id
    VolumeIdentifier getId()          const {return m_planeId; }
    //! access to the strip number
    UInt_t           getStripNumber() const {return m_strip; }
    //! access to the energy, new form
    Double_t         getEnergy()      const {return m_energy; }
    //! access to energy
    Double_t         energy()         const {return m_energy; }
    Bool_t           noise()          const {return getNoise();}
    Bool_t           getNoise()       const {return m_noise;}
    // get list of hits
    const TRefArray& getHits()        const {return m_hits;}

    //! add a hit to the list if it's not already there
    void addHit(McPositionHit* hit) {   
        TObject* inList = m_hits.FindObject(hit);
        if (inList == 0) {
            m_hits.Add(hit);
        }
    }
    
    //! add to the energy
    void operator+=(double de){m_energy += de;}

    //! access to alignment translation
    TVector3 getDelta()  const {return TVector3(m_deltaX, m_deltaY, 0.0); }
    
    Double_t getDeltaX() const {return m_deltaX; }
    Double_t getDeltaY() const {return m_deltaY; }
        
private:
    //! the plane that this belongs to
    VolumeIdentifier m_planeId;
    
    //! strip number
    UInt_t m_strip;
    
    //! total energy deposited
    Double_t m_energy;

    //! noise bit
    Bool_t m_noise;

    //! list of pointers to hits constributing to this strip
    TRefArray m_hits;

    //! Alignment offset of this strip
    Double_t m_deltaX;
    Double_t m_deltaY;

    ClassDef(McTkrStrip,1)  // Monte Carlo strip hits class
};

#endif
