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
#include "Event/Recon/TkrRecon/TkrRecInfo.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
/** 
* @class TkrVertex
*
* @brief Contains the information for a single found vertex
*
* @author The Tracking Software Group
*
* $Header$
*/

namespace Event { //Namespace
  
class TkrVertex : public TkrRecInfo, virtual public ContainedObject 
  {    
public:
    
    TkrVertex( int layer, int tower, double energy, double quality, const Ray& testRay);
   ~TkrVertex() {}

    /// Define the TkrBase routines
    double        getQuality()                       const {return m_quality;};
    double        getEnergy(TrackEnd end = Start)    const {return m_energy;}
    int           getLayer(TrackEnd end = Start)     const {return m_firstLayer;}
    int           getTower(TrackEnd end = Start)     const {return m_itower;}
    Point         getPosition(TrackEnd end = Start)  const {return m_position;}
    Vector        getDirection(TrackEnd end = Start) const {return m_direction;}
    Ray           getRay(TrackEnd end = Start)       const {return Ray(getPosition(),getDirection());}
    TkrFitPar     getTrackPar(TrackEnd end = Start)  const {return m_vertexPar;}
    double        getTrackParZ(TrackEnd end = Start) const {return m_position.z();}
    TkrFitMatrix  getTrackCov(TrackEnd end = Start)  const {return m_vertexCov;}
    bool          empty(int numHits)                 const {return m_firstLayer >= 0;}

    // Add tracks to the list
    void addTrack(TkrFitTrack* pTrack) {m_tracks.push_back(pTrack);}
    
    // How many tracks in the vertex?
    int  getNumTracks()                {return m_tracks.size();}

    // Pointers to track info
    TkrFitTrackCol::iterator getTrackIterBegin()   {return m_tracks.begin();}
    TkrFitTrackCol::iterator getTrackIterEnd()     {return m_tracks.end();}

    /// Utilities 
    void writeOut(MsgStream& log) const; 
    
private:
    TkrFitPar    m_vertexPar;
    TkrFitMatrix m_vertexCov;
    Point        m_position;
    Vector       m_direction;
    double       m_energy;
    double       m_quality;
    int          m_firstLayer;
    int          m_itower;        
    
    TkrFitTrackCol    m_tracks;
};

//typedef for the Container
typedef ObjectVector<TkrVertex>     TkrVertexCol;

}; //Namespace

#endif
