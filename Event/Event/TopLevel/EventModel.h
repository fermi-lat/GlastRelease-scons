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

// The following is a snippet taken from the libApiExports.h file that 
// was authored by Matt Langston.
// The intent is to define for windows those classes we would like 
// to export (or import) from the EventLib dll. 
#if (defined(_WIN32) && defined(_MSC_VER))
# ifdef EVT_DLL_EXPORTS
#  undef  DLL_EXPORT_EVT
#  define DLL_EXPORT_EVT __declspec(dllexport)
# else
#  undef  DLL_EXPORT_EVT
#  define DLL_EXPORT_EVT __declspec(dllimport)
# endif
#else
// The gcc compiler (i.e. the Linux/Unix compiler) exports the Universe
// of symbols from a shared library, meaning that we can't control the
// EVT of our shared libraries. We therefore just define the Symbol
// Export Macro to expand to nothing.
# undef  DLL_EXPORT_EVT
# define DLL_EXPORT_EVT
#endif

// Define a class to hold all of our definitions
class DLL_EXPORT_EVT EventModel
{
public:
    class DLL_EXPORT_EVT MC 
    {
    public:
        MC() {}
       ~MC() {}

        static std::string Event;
        static std::string McParticleCol;
        static std::string McPositionHitCol;
        static std::string McIntegratingHitCol;
        static std::string McTrajectoryCol;
        static std::string McTkrStripCol;
        static std::string D2EntryCol;
        static std::string ExposureCol;
        static std::string McEventStructure;
        static std::string McPartToPosHitTab;
        static std::string McPartToClusTab;
        static std::string McPartToClusHitTab;
        static std::string McPartToTkrCandHitTab;
        static std::string McPartToTkrPatCandTab;
        static std::string McPartToTkrTrackHitTab;
        static std::string McPartToTkrTrackTab;
	static std::string McAcdTkrPointCol;
	static std::string McAcdTkrHitPocaCol;
	static std::string McAcdTkrAssocCol;
    };

    class DLL_EXPORT_EVT Digi
    {
    public:
        Digi() {}
       ~Digi() {}

        static std::string Event;
        static std::string AcdDigiCol;
        static std::string TkrDigiCol;
        static std::string CalDigiCol;
        static std::string CalDigiHitTab;
        static std::string TkrDigiHitTab;
        static std::string TkrClusterHitTab;
    };

    class DLL_EXPORT_EVT Overlay
    {
    public:
        Overlay() {}
       ~Overlay() {}

        static std::string Event;
        static std::string EventHeader;
        static std::string TriRowBits;
        static std::string Time;
        static std::string EventSummary;
        static std::string Gem;
        static std::string Error;
        static std::string Diagnostic;
        static std::string ObfFilterStatus;
        static std::string ObfFilterTrack;
        static std::string MetaEvent;
        static std::string Ccsds;
        static std::string AncillaryEvent;
        static std::string AncillaryEventDigi;
        static std::string AcdDigiCol;
        static std::string TkrDigiCol;
        static std::string CalDigiCol;
        static std::string CalDigiHitTab;
        static std::string TkrDigiHitTab;
        static std::string TkrClusterHitTab;
    };


    class DLL_EXPORT_EVT TkrRecon
    {
    public:
        TkrRecon() {}
       ~TkrRecon() {}

        static std::string Event;
        static std::string TkrIdClusterMMap;
        static std::string TkrIdClusterMap;
        static std::string TkrClusterCol;
        static std::string TkrTrackCol;
        static std::string TkrTrackHitCol;
        static std::string TkrVertexCol;
        static std::string TkrDiagnostics;
        static std::string TkrEventParams;
        static std::string TkrFilterParamsCol;
        static std::string TkrFilterParamsToBoxTab;
        static std::string TkrFilterParamsToPointsTab;
        static std::string TkrTruncatedPlane;
        static std::string TkrTruncationInfo;
        static std::string TkrBoundBoxCol;
        static std::string TkrBoundBoxPointsCol;
        static std::string TkrBoundBoxPointsToBoxTab;
        static std::string TkrVecPointCol;
        static std::string TkrVecPointInfo;
        static std::string TkrVecPointsLinkCol;
        static std::string TkrTrackElementsCol;
        static std::string TkrTrackElemsToLinksTab;
    };


    class DLL_EXPORT_EVT CalRecon
    {
    public:
       CalRecon() {}
      ~CalRecon() {}

       static std::string Event;
       static std::string CalXtalRecCol;
       //@@@FP 07/09/05
       // reste des premieres versions
       //            static std::string CalMIPsCol;
       //@@@FP 07/09/05
       static std::string CalMipTrackCol;
       static std::string CalClusterCol;
       static std::string CalEventEnergyCol;
       static std::string CalClusterHitTab;
       //@@@FP 07/09/05
       // idem
//       static std::string CalXtalMIPsTab;
    //@@@FP 07/09/05
           //@@@CL 06/01/06 BEGIN:
       static std::string GcrXtalCol;
       static std::string GcrReconVals;
       //@@@CL 06/01/06 END
       //@@@CL 26/27/06 BEGIN:
       static std::string GcrSelectedXtalsCol;
       //@@@CL 06/27/06 END
       static std::string GcrTrack;
       static std::string GcrSelectVals;
 
    };


    class DLL_EXPORT_EVT AcdRecon
    {
    public:
        AcdRecon() {}
       ~AcdRecon() {}

        static std::string Event;
    };

    class DLL_EXPORT_EVT AcdReconV2
    {
    public:
        AcdReconV2() {}
       ~AcdReconV2() {}

        static std::string Event;
    };

    EventModel() {}
   ~EventModel() {}

    static std::string EventHeader;
};

#endif // _H_GlastEvent_EventModel_
