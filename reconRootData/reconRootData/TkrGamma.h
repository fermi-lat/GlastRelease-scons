#ifndef TKRGAMMA_H
#define TKRGAMMA_H

#include "TObject.h"
#include "TObjArray.h"

// ===================================
//          TkrGamma
// ===================================
//  
//! Root object that represents a reconstructed gamma.
/*! A gamma contains two lists. 
  1) A list of TkrLocator objects that contain the position and angle of the 
     track at different interesting positions 
     (e.g.: initial point, exit point of the tracker).
  2.) A list of two TkrTrack objects associated with this gamma, the
      e-/e+ tracks.  */
//
//

class TkrGamma : public TObject
{


private:

    //! list of TkrLocators associated with the Gamma
    TObjArray *m_locator;   

    //! pointers to e-/e+ tracks 
    TObjArray *m_tracks; 

public:

    // constructors
    //! default constructor
    /*! client must call Create() to setup TObjArray properly*/
    TkrGamma();
    //! non-default constructor that calls Create for client!
    TkrGamma(Int_t dummy);

    //! destructor
    virtual ~TkrGamma();

    void Clean();
    void Create();

    //! provide access to the list of tracks
    // allows the caller to modify the TObjArray
    TObjArray *getTracks() { return m_tracks; };

    //! provide access to the list of hits
    // allow the caller to modify the TObjArray
    TObjArray *getLocator() { return m_locator; };

    ClassDef(TkrGamma,1)
};

#endif
