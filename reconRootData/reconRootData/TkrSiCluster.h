#ifndef ROOT_TkrSiCluster_H
#define ROOT_TkrSiCluster_H

#include "TObject.h"

/** @class TkrSiCluster
 * @brief Root Tree Object storing a single Silicon Cluster
 *   This class stores one such cluster.
 *   A cluster is defined by:
 * - A layer number
 * - center strip ID
 * - number of strips in the cluster 
 * $Header$
*/
class TkrSiCluster : public TObject
{

public:
    //! denotes type of Si plane - x or y
    typedef enum {
        X = 0,
        Y
    } TKRAxes;
 

public:
    
    // in this case there are no pointers in the class, so
    // users may utilize either constructor
    TkrSiCluster();

    virtual ~TkrSiCluster();

    //! resets all member variables to default values
    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    void initialize(UInt_t id, UInt_t layer, TKRAxes xy, UShort_t center, 
        UShort_t numStrips, Float_t pos, Float_t zpos);

    //! access the ID # associated with this TkrSiCluster
    UInt_t getId() const { return m_id; };

    //! access the layer number
    UInt_t getLayer() const { return m_layer; };
    
    //! access to the type of layer - is it an X or a Y?
    TKRAxes getXY() const { return (m_xy ? Y : X); };

    //! access the Center strip number
    UShort_t getCenterStrip() const { return m_centerStrip; };

    //! access number of strips in this cluster
    UShort_t getNumStrips() const { return m_numStrips; };

    //! access the position
    Float_t getPosition() const { return m_position; };

    //! access the z position
    Float_t getZPosition() const { return m_z; };


private :
  //! id of the cluster
  UInt_t m_id;   
  //! layer number
  UInt_t m_layer;

  //! strip position
  Float_t m_position; 

  //! z position   
  Float_t m_z;              

  //! denotes whether this is an X or Y plane
  TKRAxes m_xy;             
  
  //! center strip of the cluster
  UShort_t m_centerStrip; 
  //! number of strips in this cluster
  UShort_t m_numStrips;        
  
  ClassDef(TkrSiCluster,2) 
};
#endif
