#ifndef __TKRTRACKTAB_H
#define __TKRTRACKTAB_H 1

#include "GaudiKernel/MsgStream.h"
#include "Event/RelTable/RelTable.h"
#include "Event/RelTable/Relation.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrFitTrackBase.h"

/** 
* @class TkrTrackTab
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
    typedef Event::RelTable<Event::TkrPatCand, Event::TkrFitTrackBase>              TkrFitTrackTab;
    typedef Event::Relation<Event::TkrPatCand, Event::TkrFitTrackBase>              TkrFitTrackRel;
    typedef ObjectList< Event::Relation<Event::TkrPatCand,Event::TkrFitTrackBase> > TkrFitTrackTabList;
};

#endif
