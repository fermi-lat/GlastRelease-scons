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

        enum {VERSION = 2};

        // fields and shifts of the status word, which together make a mask
        //  

        enum { 
            fieldUSED        = 1,    // tells whether cluster is used on a track
            fieldEND         = 3,    // identifies controller, 0, 1, 2=mixed
            fieldRAWTOT      = 0xff, // raw ToT
            fieldPLANEOFFSET = 1,    // to calculate Plane number from Tray/Face (1 for LAT)
            fieldLAYEROFFSET = 1,    // to calculate Layer number from Plane (0 for LAT)
            fieldVERSION     = 0xf   // version of class, room for 15!
        };
        enum {    
            shiftUSED        =  0, 
            shiftEND         =  1,
            shiftRAWTOT      =  8,
            shiftPLANEOFFSET = 26,
            shiftLAYEROFFSET = 27,
            shiftVERSION     = 28 
        };
        enum {
            maskUSED        = fieldUSED<<shiftUSED,
            maskEND         = fieldEND<<shiftEND,
            maskRAWTOT      = fieldRAWTOT<<shiftRAWTOT,
            maskPLANEOFFSET = fieldPLANEOFFSET<<shiftPLANEOFFSET,
            maskLAYEROFFSET = fieldLAYEROFFSET<<shiftLAYEROFFSET,
            maskVERSION     = fieldVERSION<<shiftVERSION
        };
    
    TkrCluster();
	
    TkrCluster(commonRootData::TkrId tkrId, UInt_t istrip0, UInt_t istripf, 
               const TVector3& position, Double_t ToT, UInt_t status);

    virtual ~TkrCluster();
	
    /// resets all member variables to default values
    void Clear(Option_t *option="");
	
    void Print(Option_t *option="") const;
		
    inline commonRootData::TkrId getTkrId()      const { return m_tkrId;}
    inline UInt_t                getFirstStrip() const { return m_strip0; }
    inline UInt_t                getLastStrip()  const { return m_stripf; }
    inline Double_t              getToT()        const { return m_ToT; }
    //inline UInt_t                getId()         const { return m_id; }

    /// Returns the chip number calculated using the strip number
    inline Int_t getChip()     const { return (m_strip0/64); };
    /// Returns the center strip of this cluster
    inline Double_t getStrip() const { return (0.5*(m_strip0+m_stripf)); };
    /// Returns the size of this cluster
    inline Double_t getSize()  const { return (TMath::Abs(Double_t(m_stripf-m_strip0)) + 1.); };
    
    inline const TVector3& getPosition() const { return m_position; }
    /// construct plane from tray/face
    inline UInt_t getStatusWord() const { return m_status; }
    inline int getPlane() const {
        return m_tkrId.getTray()*2 + m_tkrId.getBotTop() - getPlaneOffset(); 
    }
    // construct layer from Plane
    inline getLayer()     const { return (getPlane() + getLayerOffset())/2; }
    
    /// returns true if the cluster has been flagged
        bool hitFlagged()     const { 
            return ((maskVERSION&m_status)!=0 ? ((m_status&maskUSED)>0) : m_status!=0);}
	
private:
        inline int getPlaneOffset() const { 
            return ((maskVERSION & m_status)!=0 ? ((m_status&maskPLANEOFFSET)>>shiftPLANEOFFSET) : 1); }
        inline int getLayerOffset() const { 
            return ((maskVERSION & m_status)!=0 ?((m_status&maskLAYEROFFSET)>>shiftLAYEROFFSET) : 0); }
    
    /// volume id
    commonRootData::TkrId m_tkrId;  
    /// initial strip address of the cluster
    UInt_t   m_strip0;
    /// final strip address of the cluster
    UInt_t   m_stripf;
    /// space position of the cluster
    TVector3 m_position;
    /// ToT value of the cluster
    Double_t m_ToT;    
    /// Everything  below here is a candidate for removal
    /// flag of the cluster, used during pattern recognition
    UInt_t   m_status;
    //TO BE REMOVED
    UInt_t   m_id;   

	ClassDef(TkrCluster,2) 
};

#endif
