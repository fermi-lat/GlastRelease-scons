//----------------------------------------
//
//      Contains fit track information for a given plane
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections
//      W.B. Atwood, SCIPP/UCSC, Nov. 2001   
//      
//----------------------------------------

#ifndef _TkrFitPlane_H
#define _TkrFitPlane_H 1

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitHit.h"

/** 
* @class TkrFitPlane
*
* @brief Contains the data members which specify hit information 
*        for a single plane in a TKR fit track, and access methods
*
* Adapted from original KalPlane version of Bill Atwood
*
* @author Bill Atwood
*
* $Header$
*/

namespace Event { // Namespace

class TkrFitPlane
{   // Class to link together the various hit types associated with 
    // each measuring plane.  
    // Note: X and Y planes are recored individually

public:
    typedef TkrCluster::view AXIS;
    
    TkrFitPlane() 
        : m_IDHit (0), m_IDTower(0),  m_IDPlane(0), m_projection(TkrCluster::XY), 
		  m_projPlus(TkrCluster::XY), m_zplane(0), m_eneplane(0), m_radLen(0),  
		  m_activeDist(0)
    {}

    TkrFitPlane(unsigned id, int tower, int kplane, double ene, double z, const TkrFitHit& hit, AXIS prj)
        : m_IDHit(id), m_IDTower(tower), m_IDPlane(kplane), m_projection(prj), m_projPlus(TkrCluster::XY), 
		  m_zplane(z), m_eneplane(ene), m_radLen(0), m_activeDist(0)   
    { 
	    setHit(hit);
    }
    TkrFitPlane(unsigned id, int tower, int kplane, double ene, double z, AXIS prj)
        : m_IDHit(id), m_IDTower(tower), m_IDPlane(kplane), m_projPlus(TkrCluster::XY), m_eneplane(ene),
          m_radLen(0), m_activeDist(0)
    {
        m_zplane = z;
        m_projection = prj;
    }

    void     initializeInfo(unsigned int hit, unsigned int tower, unsigned int plane, AXIS proj,
                            AXIS nextProj, double z, double energy, double radLen, double activeDist);

    void     initializeHits(const TkrFitHit& meas, const TkrFitHit& pred, const TkrFitHit& fit,
                            const TkrFitHit& smooth, const TkrFitMatrix& material);


    // Access Information
    inline unsigned getIDHit()                           const {return m_IDHit;}
    inline int      getIDTower()                         const {return m_IDTower;}
    inline int      getIDPlane()                         const {return m_IDPlane;}
    inline double   getZPlane()                          const {return m_zplane;}
    inline double   getEnergy()                          const {return m_eneplane;}
    inline AXIS     getProjection()                      const {return m_projection;}
    inline AXIS     getNextProj()                        const {return m_projPlus;}
    inline double   getRadLen()                          const {return m_radLen;}
    TkrFitHit       getHit(TkrFitHit::TYPE type)         const;
    Point           getPoint(TkrFitHit::TYPE type)       const;
    double          getDeltaChiSq(TkrFitHit::TYPE type)  const;
    double          getDeltaChiEne(TkrFitHit::TYPE type) const; 
    double          getSigma(TkrFitHit::TYPE type)       const;
    TkrFitMatrix    getQmaterial()                       const {return m_Qmaterial;} 
    inline double   getActiveDist()                      const {return m_activeDist;}

    /// Classes allowed access for filling the information
    friend class KalFitTrack;
    friend class KalFitter;
    friend class TkrKalFitTrack;
    friend class KalmanFilter;

protected:
    // These methods are only accessible to friend classes for filling 
    // Set methods for filling information
    // Adding Hits, energy and other variables
    void            setHit(const TkrFitHit& hit);
    inline void     setEnergy(double const e)                  {m_eneplane = e;}
    inline void     setIDHit(unsigned id)                      {m_IDHit = id; 
                                                                m_IDTower = (int) id/1000000;}
    inline void     setIDPlane(int id)                         {m_IDPlane = id;}
    inline void     setZPlane(double z)                        {m_zplane = z;}
    inline void     setRadLen(double rl)                       {m_radLen = rl;}

    inline void     setActiveDist(double d)                    {m_activeDist = d;}
    inline void     setNextProj(AXIS nextProj)                 {m_projPlus = nextProj;}
    inline void     setQmaterial(const TkrFitMatrix& Q)        {m_Qmaterial = Q;}

    void setDeltaEne(double ene);
    void clean();   // clean the PRED - FIT - SMOOTH values but not the MEAS
    void clear();   // clean everything
    void removeHit();
   
private:
    
    unsigned     m_IDHit;      // SiCluster Index - code = tower+layer+hit
    int          m_IDTower;    // Tower number
    int          m_IDPlane;    // Plane number
    AXIS         m_projection; // X or Y measuring SSD plane
    AXIS         m_projPlus;   // X or Y measuring SSD plane+1
    double       m_zplane;     // Z location of plane
    double       m_eneplane;   // Energy of track at this plane
    double       m_radLen;     // Rad. Len. associated with this palne
    double       m_activeDist; // The distance inside (positive) hit SSD (neg. if outside)
        
    TkrFitHit    m_hitmeas;
    TkrFitHit    m_hitpred;
    TkrFitHit    m_hitfit;
    TkrFitHit    m_hitsmooth;
    TkrFitMatrix m_Qmaterial;  // covarianve matrix of the material effects 
};

//Following typedefs for containing fit track objects
typedef std::vector<TkrFitPlane>                  TkrFitPlaneCol;
typedef std::vector<TkrFitPlane>::iterator        TkrFitPlaneColPtr;
typedef std::vector<TkrFitPlane>::const_iterator  TkrFitPlaneConPtr;

bool operator<(const TkrFitPlane&, const TkrFitPlane&); 
bool operator==(const TkrFitPlane&, const TkrFitPlane&); 

}; //Namespace

#endif
