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
        };

        namespace Irf {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string IrfTkrHits;
            _EXTERN_ std::string IrfCalHits;
            _EXTERN_ std::string IrfAcdHits;
        }

        namespace Raw {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string TdCsIDatas;
        };


    };
#undef _EXTERN_

#endif // GLASTEVENT_EVENTMODEL_H
