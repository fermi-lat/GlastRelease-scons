#ifndef __TkrPatCand_H
#define __TkrPatCand_H
/** 
* @class TkrPatCand
*
* @brief Contains the candidate track information needed for the track fit
*        Note that everything needed to seed the fitter (in find hits mode)
*        can be obtained by the methods defined in TkrRecInfo.
*
* @author The Tracking Software Group
*
* $Header$
*/
#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "Event/Recon/TkrRecon/TkrRecInfo.h"
#include "Event/Recon/TkrRecon/TkrPatCandHit.h"

extern const CLID& CLID_TkrPatCand;

namespace Event { //Namespace

class TkrPatCand: public TkrRecInfo, virtual public ContainedObject
{    
public:
    
    TkrPatCand(int layer, int tower, double energy, double energyErr, double quality, const Ray& testRay);
    TkrPatCand(int layer, int tower, double energy, double quality, const Ray& testRay) 
                    {TkrPatCand(layer,tower,energy,1.,quality,testRay);}
    ~TkrPatCand();

   //! Retrieve pointer to class defininition structure
   virtual const CLID& clID() const   { return TkrPatCand::classID(); }
   static const CLID& classID()       { return CLID_TkrPatCand; }

    /// Define the TkrRecInfo routines

    double        getQuality()                       const {return m_quality;};
    double        getEnergy(TrackEnd end = Start)    const ;
    double        getEnergyErr(TrackEnd end = Start) const ;
    int           getLayer(TrackEnd end = Start)     const ;
    int           getTower(TrackEnd end = Start)     const ;
    Point         getPosition(TrackEnd end = Start)  const ;
    Vector        getDirection(TrackEnd end = Start) const ;
    Ray           getRay(TrackEnd end = Start)       const ;
    TkrFitPar     getTrackPar(TrackEnd end = Start)  const ;
    double        getTrackParZ(TrackEnd end = Start)   const ;
    TkrFitMatrix  getTrackCov(TrackEnd end = Start)    const ; 
    bool          empty(int numHits)                const ;

    //Provide a method to writeout the contents of the class
    void writeOut(MsgStream& log) const; 

    //Method to store hits into vector
    void addCandHit(TkrCluster*   pCluster);
    void addCandHit(TkrPatCandHit candHit);

    //Methods to set member variables if need be
    void setEnergy(double energy)       {m_energy = energy;}
    void setEnergyErr(double energyErr) {m_energyErr = energyErr;}

    //Access to some information regarding hits (if they exist)
    int              lastLayer();

    //Provide access to the vector of hit candidates
    int              numPatCandHits()    {return m_hits.size();}
    TkrPatCandHit*   getCandHit(int idx) {return m_hits[idx];}
    TkrPatCandHit*   getFoLPlane(TrackEnd end = Start);

    //Provide access to a vector iterator (do we want this?)
    CandHitVectorPtr getHitIterBegin()   {return m_hits.begin();}

    //Provide an stl like interface to the hits
    CandHitVectorPtr begin()             {return m_hits.begin();}
    CandHitVectorPtr end()               {return m_hits.end();}
    int              size()              {return m_hits.size();}

    //Provide mechanism to sort hits in the list
    void             sortHits();
           
private:
    //For sorting the hits
    ////friend bool      sortHits(const TkrPatCandHit o1, const TkrPatCandHit o2) {return o1 < o2;}

    CandHitVector m_hits;

    //Information mandated by TkrRecInfo
    Point         m_position;
    Vector        m_direction;
    double        m_energy;
    double        m_energyErr;
    double        m_quality;
    int           m_firstLayer;
    int           m_itower;    
};

//typedef for the Container
typedef ObjectVector<TkrPatCand>      TkrPatCandCol;
typedef TkrPatCandCol::const_iterator TkrPatCandColConstPtr;
typedef TkrPatCandCol::iterator       TkrPatCandColPtr;

}; //Namespace

#endif
