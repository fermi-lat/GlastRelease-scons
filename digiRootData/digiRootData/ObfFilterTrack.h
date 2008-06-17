/** @file ObfFilterTrack.h
* @author Tracy Usher
*
* $Header$

*/
#ifndef ObfFilterTrackRoot_H
#define ObfFilterTrackRoot_H 1

#include "TObject.h"
#include "TObjArray.h"

/** 
* @class ObfFilterTrack
*
* @brief This TDS object contains the "best" track information from 
*        the onboard filter tracking.
*
*/
class ObfFilterTrack : public TObject
{
public:

    // Standard constructor
    ObfFilterTrack() : TObject(), m_nXhits(0), m_nYhits(0), m_xInt(0), m_yInt(0), m_z(0), m_slpXZ(0), m_slpYZ(0) {}

    // Constructor with known parameters
    ObfFilterTrack(Int_t xHits, Int_t yHits, Float_t xInt, Float_t yInt, Float_t z, Float_t slpXZ, Float_t slpYZ) : TObject(),
        m_nXhits(xHits), m_nYhits(yHits), m_xInt(xInt), m_yInt(yInt), m_z(z), m_slpXZ(slpXZ), m_slpYZ(slpYZ) {}

    // Standard destructor
    virtual ~ObfFilterTrack() {}

    // A special "fill" method for use in Root IO
    void initialize(Int_t   xHits = 0, 
                    Int_t   yHits = 0, 
                    Float_t xInt  = 0., 
                    Float_t yInt  = 0., 
                    Float_t z     = 0., 
                    Float_t slpXZ = 0., 
                    Float_t slpYZ = 0.);

    // Get the track paramters 
    const Int_t   get_nXhits() const {return m_nXhits;}
    const Int_t   get_nYhits() const {return m_nYhits;}
    const Float_t get_xInt()   const {return m_xInt;}
    const Float_t get_yInt()   const {return m_yInt;}
    const Float_t get_z()      const {return m_z;}
    const Float_t get_slpXZ()  const {return m_slpXZ;}
    const Float_t get_slpYZ()  const {return m_slpYZ;}

    /// clear lists, free pointers, etc., after read from / write to file
    void   Clear(Option_t *option ="");
    void   Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const ObfFilterTrack &, const std::string & name = "" ) const ; // for tests
    void   Print(Option_t *option="") const;

private:
    Int_t   m_nXhits;  // Number of hits on the X projection
    Int_t   m_nYhits;  // Number of hits on the Y projection
    Float_t m_xInt;    // Intercept of the X projection
    Float_t m_yInt;    // Intercept of the Y projection
    Float_t m_z;       // Z coordinate of intercept
    Float_t m_slpXZ;   // slope X-Z plane
    Float_t m_slpYZ;   // slope Y-Z plane
    
    ClassDef(ObfFilterTrack,1) // Onboard Filter Track class
};


#endif // ObfFilterTrack_H
