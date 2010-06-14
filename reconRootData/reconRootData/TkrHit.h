#ifndef TKRHIT_H
#define TKRHIT_H

#include "TObject.h"

class TkrSiCluster;
class TkrLocator;

//! Root object concerning a silicon cluster that is in track.
/*! This object is only referenced via a TkrTrack
/*! Contains a pointer to a TkrSiCluster that is associated with this TkrHit.  
    Also contains stats pertaining to this TkrHit within a TkrTrack. 
*/
    

class TkrHit : public TObject 
{

private :
         
  Float_t m_residual;       // residual of the silicon cluster
  Float_t m_chi2;           // chi**2 

  TkrSiCluster *m_cluster;  // pointer to the TkrSiCluster 
                            //  corresponding
                            // to this TkrHit - part of the
                            // TkrRecon list of TkrSiClusters

  TkrLocator *m_locator;    // pointer to a TkrLocator

public :

    //! constructor
    TkrHit();

    //! destructor
    virtual ~TkrHit();

    void Clean();
    void Create();

    //! provide access to the residual
    Float_t getResidual() { return m_residual; };
    void setResidual(Float_t r) { m_residual = r; };

    //! access Chi2
    Float_t getChi2() { return m_chi2; };
    void setChi2(Float_t chi) { m_chi2 = chi; };  

    //! acces the TkrSiCluster pointer
    TkrSiCluster* getSiCluster() { return m_cluster; };
    //! set the TkrSiCluster pointer
    inline void setSiCluster(TkrSiCluster *si) { m_cluster = si; };

    //! access TkrLocator pointer
    TkrLocator* getLocator() { return m_locator; };
    //! set the TkrLocator pointer
    inline void setLocator(TkrLocator *l) { m_locator = l; };

    ClassDef(TkrHit,1) 
};
#endif
