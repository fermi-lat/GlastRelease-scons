#ifndef TkrPatCandCol_h
#define TkrPatCandCol_h
/** 
* @class TkrPatCandCol
*
* @brief TDS Container class to hold the list of candidate tracks
*
* @author The Tracking Software Group
*
* $Header$
*/

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_TkrPatCandCol;

namespace TkrRecon { //Namespace

class TkrPatCandCol : public DataObject
{
public:
    TkrPatCandCol();
   ~TkrPatCandCol();

    //Provide ability to output some information about self...
    void writeOut(MsgStream& log) const;

	//! GAUDI members to be use by the converters
	static const CLID&  classID()           {return CLID_TkrPatCandCol;}
	virtual const CLID& clID()        const {return classID();}

    //How many track candidates are there?
    int                 getNumCands() const {return m_candTracks.size();}

    //Access to tracks through an iterator
    CandTrkVectorPtr    getTrackPtr()       {return m_candTracks.begin();}

    //Access to tracks by index
    TkrPatCand*         getTrack(int idx)   {return m_candTracks[idx];}

    //Add tracks to the list
    void                addTrack(TkrPatCand* candTrack);

private:
    CandTrkVector m_candTracks;
};

}; //Namespace

#endif