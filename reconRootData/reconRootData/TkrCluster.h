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
    
    TkrCluster();
	
    TkrCluster(commonRootData::TkrId tkrId, UInt_t istrip0, UInt_t istripf, 
               const TVector3& position, Double_t ToT, UInt_t flag, UInt_t id);

    virtual ~TkrCluster();
	
    /// resets all member variables to default values
    void Clear(Option_t *option="");
	
    void Print(Option_t *option="") const;
		
    inline commonRootData::TkrId getTkrId()      const { return m_tkrId;}
    inline UInt_t                getFirstStrip() const { return m_strip0; }
    inline UInt_t                getLastStrip()  const { return m_stripf; }
    inline Double_t              getToT()        const { return m_ToT; }
    inline UInt_t                getId()         const { return m_id; }

    /// Returns the chip number calculated using the strip number
    inline Int_t getChip()     const { return (m_strip0/64); };
    /// Returns the center strip of this cluster
    inline Double_t getStrip() const { return (0.5*(m_strip0+m_stripf)); };
    /// Returns the size of this cluster
    inline Double_t getSize()  const { return (TMath::Abs(Double_t(m_stripf-m_strip0)) + 1.); };
    
    inline const TVector3& getPosition() const { return m_position; }
    
    /// returns true if the cluster has been flagged
    inline bool hitFlagged() const { return (m_flag!=0); }
	
private:
    
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
    UInt_t   m_flag;
    //TO BE REMOVED
    UInt_t   m_id;   

	ClassDef(TkrCluster,2) 
};

#endif
