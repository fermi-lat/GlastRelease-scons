/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include "Event/Recon/TkrRecon/TkrFitTrackCol.h"

using namespace TkrRecon;

TkrFitTrackCol::~TkrFitTrackCol()
{
    int          numTracks = getNumTracks();
    TkrFitColPtr trks      = getTrackPtr();

    while(numTracks--) delete *trks++;
}


void TkrFitTrackCol::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrFitTrackCol::writeOut --- " << endreq;
    log << MSG::DEBUG << "     Number of tracks: " << getNumTracks() << endreq;

    int trackNo = 0;

    while(trackNo < getNumTracks())
    {
        log << MSG::DEBUG << "     Track Number: " << trackNo << endreq;
        m_Tracks[trackNo++]->writeOut(log);
    }

    return;
}

void TkrFitTrackCol::addTrack(TkrFitTrack* Track)
{
    m_Tracks.push_back(Track);

    return;
}