/*
        Code to implement the PatRecTracks class
        Tracy Usher Nov 27, 2000
*/

#include "Event/Recon/TkrRecon/TkrPatCandCol.h"

using namespace TkrRecon;

TkrPatCandCol::TkrPatCandCol()
{
    m_candTracks.clear();

        return;
}

TkrPatCandCol::~TkrPatCandCol()
{
    int              numCands = getNumCands();
    CandTrkVectorPtr cands    = getTrackPtr();

    while(numCands--) delete *cands++;

        return;
}


void TkrPatCandCol::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrPatCandCol::writeOut --- " << endreq;
    log << MSG::DEBUG << "     Number of candidate tracks: " << getNumCands() << endreq;

    int trackNo = 0;

    while(trackNo < getNumCands())
    {
        log << MSG::DEBUG << "     Track Number: " << trackNo << endreq;
        m_candTracks[trackNo++]->writeOut(log);
    }

    return;
}

void TkrPatCandCol::addTrack(TkrPatCand* candTrack)
{
    m_candTracks.push_back(candTrack);

    return;
}