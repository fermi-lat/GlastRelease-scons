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

// Include files
#include <string>

#if defined(_GlastEvent_EventModel_CPP_)
#define  _EXTERN_ 
#else
#define  _EXTERN_ extern
#endif

/*  The following names will be corrected in next releases
 *    McParticle    becomes   McParticles
 *    McVertex      becomes   McVertices
 *    AcdDigi       becomes   AcdDigis
 */  
    namespace EventModel {
        _EXTERN_ std::string   Event;

        namespace MC {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string McVertices;
            _EXTERN_ std::string McParticles;
            _EXTERN_ std::string McPositionHits;
            _EXTERN_ std::string McIntegratingHits;
        }

        namespace Digi {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string AcdDigis;
        }

        namespace Irf {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string IrfTkrHits;
            _EXTERN_ std::string IrfCalHits;
            _EXTERN_ std::string IrfAcdHits;
        }

        namespace Raw {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string TdSiDatas;
            _EXTERN_ std::string TdCsIDatas;
        }

        namespace TkrRecon {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string SiLayers;
            _EXTERN_ std::string SiClusters;
            _EXTERN_ std::string SiRecObjs;
        }
    }

#undef _EXTERN_
#endif // GLASTEVENT_EVENTMODEL_H
