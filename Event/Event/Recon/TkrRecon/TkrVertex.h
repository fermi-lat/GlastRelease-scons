//
//-----------------------------------------------------------------
//
//  TkrVertex
//
//  Class definition for a found vertex 
//  ** Test Version **
//
//  Adapted from TkrFitTrack.h
//  Tracy Usher March 1, 2002
//
//-----------------------------------------------------------------
//
#ifndef __TkrVertex_H
#define __TkrVertex_H 1

#include <vector>
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/IInterface.h"
#include "Event/Recon/TkrRecon/TkrRecInfo.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
/** 
* @class TkrVertex
*
* @brief Contains the information for a single found vertex
*
* @author The Tracking Software Group
*
* $Header$
*/

static const CLID& CLID_TkrVertex = InterfaceID("TkrVertex", 1, 0);

namespace Event { //Namespace
  
class TkrVertex : public TkrRecInfo, virtual public ContainedObject 
  {    
public:
    
    TkrVertex() {}
    TkrVertex( int layer, int tower, double energy, double quality, const Ray& testRay);
   ~TkrVertex() {}

   //! Retrieve pointer to class defininition structure
   virtual const CLID& clID() const   { return TkrVertex::classID(); }
   static const CLID& classID()       { return CLID_TkrVertex; }


    /// Define the TkrBase routines
    double        getQuality() const;                    
    double        getEnergy(TrackEnd end = Start)    const;   
    int           getLayer(TrackEnd end = Start)     const; 
    int           getTower(TrackEnd end = Start)     const; 
    Point         getPosition(TrackEnd end = Start)  const; 
    Vector        getDirection(TrackEnd end = Start) const; 
    Ray           getRay(TrackEnd end = Start)       const; 
    TkrFitPar     getTrackPar(TrackEnd end = Start)  const; 
    double        getTrackParZ(TrackEnd end = Start) const; 
    TkrFitMatrix  getTrackCov(TrackEnd end = Start)  const; 
    bool          empty(int numHits)                 const; 

    // Add tracks to the list
    void addTrack(TkrTrack* pTrack) {m_tracks.push_back(pTrack);}
    
    // How many tracks in the vertex?
    int  getNumTracks() const {return m_tracks.size();}

    // Pointers to track info
    SmartRefVector<TkrTrack>::const_iterator getTrackIterBegin() const {return m_tracks.begin();}
    SmartRefVector<TkrTrack>::const_iterator getTrackIterEnd()   const {return m_tracks.end();}

    /// Utilities 
    void writeOut(MsgStream& log) const; 
    
private:
    TkrFitPar      m_vertexPar;
    TkrFitMatrix   m_vertexCov;
    Point          m_position;
    Vector         m_direction;
    double         m_energy;
    double         m_quality;
    int            m_firstLayer;
    int            m_itower; 
    
    SmartRefVector<TkrTrack> m_tracks;
};

//typedef for the Container
typedef ObjectVector<TkrVertex>      TkrVertexCol;
typedef TkrVertexCol::const_iterator TkrVertexConPtr;
typedef TkrVertexCol::iterator       TkrVertexColPtr;

}; //Namespace

#endif
