#ifndef __TKRVERTEXTAB_H
#define __TKRVERTEXTAB_H 1

#include "GaudiKernel/MsgStream.h"
#include "Event/RelTable/RelTable.h"
#include "Event/RelTable/Relation.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrFitTrackBase.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

/** 
* @class TkrVertexTab
*
* @brief TkrRecon header file for relational table typedefs
*
*
* @author(s) The Tracking Software Group
*
* $Header$
*/

namespace Event
{
    typedef Event::RelTable<Event::TkrVertex, Event::TkrFitTrackBase>              TkrVertexTab;
    typedef Event::Relation<Event::TkrVertex, Event::TkrFitTrackBase>              TkrVertexRel;
    typedef ObjectList< Event::Relation<Event::TkrVertex,Event::TkrFitTrackBase> > TkrVertexTabList;

    typedef Event::RelTable<Event::TkrVertex, Event::TkrTrack>                     TkrVertexTrackTab;
    typedef Event::Relation<Event::TkrVertex, Event::TkrTrack>                     TkrVertexTrackRel;
    typedef ObjectList< Event::Relation<Event::TkrVertex,Event::TkrTrack> >        TkrVertexTrackTabList;
};

#endif
