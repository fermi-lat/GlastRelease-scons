// $Id$

#ifndef _H_GlastEvent_EventModel_
#define _H_GlastEvent_EventModel_

/* Definition of the event structure in the Transient Data Store.
 *  
 * Only two levels in the logical path are foreseen at present, 
 *   /event/<namespace>/<leave>  e.g. /Event/MC/McVertices
 * 
 * Convention:
 *  If the <leave> object is a
 *  DataObject    use name of corresponding class
 *  Container     use name of ContainedObject class in plural
 *                or append 'Vec' to the name, e.g. use
 *                McVertices or McVertexVec
 *                
 *
 * @author : adapted from LHCb EventModel
 * @author   I. Gable
 * @author   S. Gillespie
 * @author   T. H.-Kozanecka
 */ 

#include <string>

#if defined(_Event_EventModel_CPP_)
#define  _EXTERN_ 
#else
#define  _EXTERN_ extern
#endif

    namespace EventModel {
        _EXTERN_ std::string   EventHeader;

        namespace MC {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string McParticleCol;
            _EXTERN_ std::string McPositionHitCol;
            _EXTERN_ std::string McIntegratingHitCol;
            _EXTERN_ std::string McTkrStripCol;
        }

        namespace Digi {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string AcdDigiCol;
            _EXTERN_ std::string TkrDigiCol;
            _EXTERN_ std::string CalDigiCol;
        }

        namespace TkrRecon {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string SiLayers;
            _EXTERN_ std::string TkrClusterCol;
            _EXTERN_ std::string TkrPatCandCol;
            _EXTERN_ std::string SiRecObjs;
            _EXTERN_ std::string TkrFitTrackCol;
            _EXTERN_ std::string TkrVertexCol;
        }

        namespace CalRecon {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string CalXtalRecCol;
            _EXTERN_ std::string CalClusterCol;
		}

        namespace AcdRecon {
            _EXTERN_ std::string Event;
        }
    }

#undef _EXTERN_
#endif // GLASTEVENT_EVENTMODEL_H
