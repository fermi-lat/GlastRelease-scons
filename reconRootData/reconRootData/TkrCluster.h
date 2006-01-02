#ifndef ROOT_TkrCluster_H
#define ROOT_TkrCluster_H

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

#include "commonRootData/idents/TkrId.h"

/** @class TkrCluster
* @brief Root Tree Object storing a single Silicon Cluster
* 
*   This class stores one such cluster.
*   A cluster contains
* - A Tracker Id (TkrId) which gives tower(s), tray, view, etc.
* - First strip number associated with the cluster
* - Last strip number associated with the cluster
* - ToT value associated with the cluster
* - Position that denotes the 3D position of the cluster
*
* $Header$
*/
class TkrCluster : public TObject
{

public:
    // once we have an official release, the version number can be used
    //  to allow backward compatibility

    enum {VERSION = 3};

    // fields and shifts of the status word, which together make a mask
    //  

    enum { 
        fieldUSED        = 1,    // tells whether cluster is used on a track
        fieldEND         = 3,    // identifies controller, 0, 1, 2=mixed
        fieldPLANEOFFSET = 1,    // to calculate Plane number from Tray/Face (1 for LAT)
        fieldLAYEROFFSET = 1,    // to calculate Layer number from Plane (0 for LAT)
   };
    enum {    
        shiftUSED        =  0, 
        shiftEND         =  1,
        shiftPLANEOFFSET = 29,
        shiftLAYEROFFSET = 30,
    };
    enum {
        maskUSED        = fieldUSED<<shiftUSED,
        maskEND         = fieldEND<<shiftEND,
        maskPLANEOFFSET = fieldPLANEOFFSET<<shiftPLANEOFFSET,
        maskLAYEROFFSET = fieldLAYEROFFSET<<shiftLAYEROFFSET,
    };

    TkrCluster();

    TkrCluster(commonRootData::TkrId tkrId, UInt_t istrip0, UInt_t istripf, 
        const TVector3& position, UInt_t rawToT, Float_t ToT, UInt_t status, UInt_t nBad);

    virtual ~TkrCluster();

    /// resets all member variables to default values
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrCluster &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    inline commonRootData::TkrId getTkrId()      const { return m_tkrId;}
    inline UInt_t                getFirstStrip() const { return m_strip0; }
    inline UInt_t                getLastStrip()  const { return m_stripf; }
    inline Float_t               getToT()        const { return getRawToT(); }
    inline UInt_t                getNBad()       const { return m_nBad; }


    /// Returns the center strip of this cluster
    inline Float_t getStrip() const { return (0.5*(m_strip0+m_stripf)); };
    /// Returns the size of this cluster
    inline Float_t getSize()  const { return (TMath::Abs(Double_t(m_stripf-m_strip0)) + 1.); };

    inline const TVector3& getPosition() const { return m_position; }
    /// construct plane from tray/face
    inline UInt_t getStatusWord() const { return m_status; }
    inline int getPlane() const {
        return m_tkrId.getTray()*2 + m_tkrId.getBotTop() - 1 + getPlaneOffset(); 
    }
    // construct layer from Plane
    inline int getLayer() const { return (getPlane() + getLayerOffset())/2; }

    /// returns true if the cluster has been flagged
    bool hitFlagged()     const { return ((m_status&maskUSED)>0);}

    /// retrieves raw ToT (will be raw or corrected depending on the version)
    inline Float_t getRawToT() const { return m_rawToT; }

    /// retrieve corrected ToT (zero for old-style records)
    inline Float_t getMips() const { return   m_ToT ; }
    /// retrieve end
    inline UInt_t getEnd() const { return ((m_status&maskEND)>>shiftEND);}

private:
    inline UInt_t getPlaneOffset() const { 
        return ((m_status&maskPLANEOFFSET)>>shiftPLANEOFFSET); }
    inline UInt_t getLayerOffset() const { 
        return ((m_status&maskLAYEROFFSET)>>shiftLAYEROFFSET); }

    /// volume id
    commonRootData::TkrId m_tkrId;  
    /// initial strip address of the cluster
    UInt_t   m_strip0;
    /// final strip address of the cluster
    UInt_t   m_stripf;
    /// number of bad strips in this cluster
    UInt_t   m_nBad;
    /// space position of the cluster
    TVector3 m_position;
    /// Raw ToT
    UInt_t m_rawToT;
    /// ToT value of the cluster
    Float_t m_ToT;    
    /// flag of the cluster, used during pattern recognition
    UInt_t   m_status;

    ClassDef(TkrCluster,2) 
};

#endif
