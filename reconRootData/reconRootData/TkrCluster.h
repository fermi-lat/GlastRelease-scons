#ifndef ROOT_TkrCluster_H
#define ROOT_TkrCluster_H

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

/** @class TkrCluster
* @brief Root Tree Object storing a single Silicon Cluster
* 
*   This class stores one such cluster.
*   A cluster contains
* - Id number
* - Tower Id number
* - Plane number
* - First strip number associated with the cluster
* - Last strip number associated with the cluster
* - View denoting if the cluster measures X or Y
* - ToT value associated with the cluster
* - Position that denotes the 3D position of the cluster
*
* $Header$
*/
class TkrCluster : public TObject
{
	
public:
    /// enumeration of the view of the cluster
    enum view 
    {
        X=0, /**< cluster measures X */ 
        Y=1  /**< cluster measures Y */ 
    };

	
	
public:
    
    TkrCluster();
	
    TkrCluster(UInt_t id, UInt_t iplane, view v, 
        UInt_t istrip0, UInt_t istripf, const TVector3& position, 
		Double_t ToT, UInt_t flag, UInt_t tower = 0);

    virtual ~TkrCluster();
	
    /// resets all member variables to default values
    void Clear(Option_t *option="");
	
    void Print(Option_t *option="") const;
		
	
    inline UInt_t getTower() const { return m_tower; }
    inline UInt_t getId() const { return m_id; }
    inline UInt_t getPlane() const { return m_plane; }
    inline view getView() const {return m_view;}
    inline UInt_t getFirstStrip() const { return m_strip0; }
    inline UInt_t getLastStrip() const { return m_stripf; }
    inline UInt_t getToT() const { return m_ToT; }

    /// Returns the chip number calculated using the strip number
    inline Int_t getChip() const { return (m_strip0/64); };
    /// Returns the center strip of this cluster
    inline Double_t getStrip() const { return (0.5*(m_strip0+m_stripf)); };
    /// Returns the size of this cluster
    inline Double_t getSize() const { return (TMath::Abs(Double_t(m_stripf-m_strip0)) + 1.); };
    
    inline const TVector3& getPosition() const { return m_position; }
    
    /// returns true if the cluster has been flagged
    inline bool hitFlagged() const { return (m_flag!=0); }
    
	
private:
    /// ToT value of the cluster
    Double_t m_ToT;

	/// id of the cluster, sequential in order of construction
    UInt_t m_id;   

    UInt_t m_tower;

    /// plane id: this is the bi-layer number, but ordered backwards (0 at front)
    UInt_t m_plane;

    /// flag of the cluster, used during pattern recognition
    UInt_t m_flag;

    /// initial strip address of the cluster
    UInt_t m_strip0;
    /// final strip address of the cluster
    UInt_t m_stripf;

	/// view (X or Y)
    TkrCluster::view m_view;
    
    /// space position of the cluster
    TVector3 m_position;
	
	ClassDef(TkrCluster,1) 
};

#endif
