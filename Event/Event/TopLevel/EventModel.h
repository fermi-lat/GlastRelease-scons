// $Id$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//
#ifndef _H_GlastEvent_EventModel_
#define _H_GlastEvent_EventModel_

// Include files
#include <string>

#if defined(_GlastEvent_EventModel_CPP_)
#define  _EXTERN_ 
#else
#define  _EXTERN_ extern
#endif

    
    namespace EventModel {
        _EXTERN_ std::string   Event;

        namespace MC {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string MCTKRHits;
            _EXTERN_ std::string MCCalorimeterHits;
            _EXTERN_ std::string MCACDHits;
        };

    };
#undef _EXTERN_

#endif // GLASTEVENT_EVENTMODEL_H
