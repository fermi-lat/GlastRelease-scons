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
       static std::string CalEventEnergy;
       static std::string CalClusterHitTab;
       //@@@FP 07/09/05
       // idem
//       static std::string CalXtalMIPsTab;
    //@@@FP 07/09/05
    };


    class DLL_EXPORT_EVT AcdRecon
    {
    public:
        AcdRecon() {}
       ~AcdRecon() {}

        static std::string Event;
    };

    EventModel() {}
   ~EventModel() {}

    static std::string EventHeader;
};

#endif // _H_GlastEvent_EventModel_
