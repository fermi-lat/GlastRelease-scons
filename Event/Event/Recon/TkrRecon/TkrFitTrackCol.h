
#ifndef TkrFitTrackCol_H
#define TkrFitTrackCol_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GlastEvent/Recon/TkrRecon/TkrFitTrack.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_TkrTracks;

/** 
* @class TkrFitTrackCol
*
* @brief Gaudi TDS container class for keeping track of a list of TkrFitTrack objects
*
* Adapted from original version of Bill Atwood
*
* @author The Tracking Software Group
*
* $Header$
*/

namespace TkrRecon { //Namespace

class TkrFitTrackCol : public DataObject
{
public:
    TkrFitTrackCol() {m_Tracks.clear();}
   ~TkrFitTrackCol();

    //Provide ability to output some information about self...
    void writeOut(MsgStream& log) const;

	//! GAUDI members to be use by the converters
	static const CLID&  classID()           {return CLID_TkrTracks;}
	virtual const CLID& clID()        const {return classID();}

    //How many reconstructed tracks are there?
    int                 getNumTracks() const {return m_Tracks.size();}

    //Access to tracks through an iterator
    TkrFitColPtr        getTrackPtr()        {return m_Tracks.begin();}
    TkrFitColPtr        getTrackEnd()        {return m_Tracks.end();}

    //Access to tracks by index
    TkrFitTrack*        getTrack(int idx)    {return m_Tracks[idx];}

    //Add tracks to the list
    void                addTrack(TkrFitTrack* track);

private:
    TkrFitCol m_Tracks;
};

}; //Namespace

#endif