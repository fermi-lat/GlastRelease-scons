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
            _EXTERN_ std::string McVertexCol;
            _EXTERN_ std::string McParticleCol;
            _EXTERN_ std::string McPositionHitCol;
            _EXTERN_ std::string McIntegratingHitCol;
        }

        namespace Digi {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string AcdDigis;
            _EXTERN_ std::string TkrDigis;
        }

        namespace Data {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string TdSiData;
            _EXTERN_ std::string TdCsIData;
            _EXTERN_ std::string TdGlastData;
            _EXTERN_ std::string TdVetoData;
        }

        namespace TkrRecon {
            _EXTERN_ std::string Event;
            _EXTERN_ std::string SiLayers;
            _EXTERN_ std::string TkrClusters;
            _EXTERN_ std::string TkrPatCandCol;
            _EXTERN_ std::string SiRecObjs;
            _EXTERN_ std::string TkrFitTrackCol;
            _EXTERN_ std::string TkrVertexCol;
        }

        namespace AcdRecon {
            _EXTERN_ std::string Event;
        }
    }

#undef _EXTERN_
#endif // GLASTEVENT_EVENTMODEL_H
