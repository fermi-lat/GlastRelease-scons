// $Id$

#ifndef _OverlayEventModel_h
#define _OverlayEventModel_h

/* Definition of the event structure in the Transient Data Store 
 * for OverlayEvent TDS objects.
 * 
 * Convention:
 *  If the <leave> object is a
 *  DataObject    use name of corresponding class
 *  Container     use name of ContainedObject class in plural
 *                or append 'Vec' to the name, e.g. use
 *                McVertices or McVertexVec
 *                
 *
 * @author   T. Usher
 */ 
#include <string>

// The following is a snippet taken from the libApiExports.h file that 
// was authored by Matt Langston.
// The intent is to define for windows those classes we would like 
// to export (or import) from the EventLib dll. 
#if (defined(_WIN32) && defined(_MSC_VER))
# ifdef OEVT_DLL_EXPORTS
#  undef  DLL_EXPORT_OEVT
#  define DLL_EXPORT_OEVT __declspec(dllexport)
# else
#  undef  DLL_EXPORT_OEVT
#  define DLL_EXPORT_OEVT __declspec(dllimport)
# endif
#else
// The gcc compiler (i.e. the Linux/Unix compiler) exports the Universe
// of symbols from a shared library, meaning that we can't control the
// EVT of our shared libraries. We therefore just define the Symbol
// Export Macro to expand to nothing.
# undef  DLL_EXPORT_OEVT
# define DLL_EXPORT_OEVT
#endif

// Define a class to hold all of our definitions
class DLL_EXPORT_OEVT OverlayEventModel
{
public:

    class DLL_EXPORT_OEVT Overlay
    {
    public:
        Overlay() {}
       ~Overlay() {}

        static std::string EventOverlay;
        static std::string SrcOverlay;
        static std::string AcdOverlayCol;
        static std::string TkrOverlayCol;
        static std::string CalOverlayCol;
        static std::string GemOverlay;
        static std::string DiagDataOverlay;
        static std::string PtOverlay;
    };

    OverlayEventModel() {}
   ~OverlayEventModel() {}

    static std::string OverlayEventHeader;
};

#endif // _OverlayEventModel_
