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

        _EXTERN_ std::string ACDTilesName;

        namespace MC {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string MCTKRHits;
            _EXTERN_ std::string MCCalorimeterHits;
            _EXTERN_ std::string MCACDHits;
        };

        namespace Hits  {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string Glast;
            //_EXTERN_ std::string TowerName;
            //_EXTERN_ std::string SiLayersName;
            _EXTERN_ std::string ACDTilesName;
            //_EXTERN_ std::string CalorimeterName;
        
        };
    };
#undef _EXTERN_

#endif // GLASTEVENT_EVENTMODEL_H
