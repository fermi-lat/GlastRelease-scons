/**@file main.cxx
@brief main program for application to create GLAST classification trees
$Header$
*/

#include "ClassifyCal.h"
#include "ClassifyVertex.h"
#include "ClassifyCore.h"
#include "ClassifyGamma.h"


/** @page applications classify
Run one or more classifications from MeritTuple root files

classify [rootpath] [treepath] [case]

@param rootpath path to the root data
@param treepath path to tree parameters, place to put results
@param case One of: all, goodcal, vertex_thin, vertex_thick, psf_thin_vertex, 
                    psf_thick_vertex, psf_thin_track, psf_thick_track, gamma

If no args, check for the environment variables ROOTPATH, TREEPATH, and CLASSIFY_TYPE
Otherwise use development defaults.

*/
int main(int argc , char * argv[])
{
    int rc=0;
    std::string // defaults for development
        rootpath("D:\\common\\DC2\\"),
        treepath("..\\data"),
        name("all");
    if( argc==1) {
        // no args: check env vars
        const char* env = ::getenv("ROOTPATH");
        if (env!=0) rootpath = std::string(env);
        env = ::getenv("TREEPATH");
        if( env!=0) treepath = std::string(env);
        env = ::getenv("CLASSIFY_TYPE");
        if( env!=0) name = std::string(env);
    }else {
        if( argc>1) rootpath=argv[1];
        if( argc>2) treepath=argv[2];
        if( argc>3) name   = argv[3];
    }
    bool all(name=="all"); 
    bool psf(name=="psf");

    GlastClassify::setPaths(rootpath, treepath);

    std::cout << "classifier invoked with root, tree paths: "<< rootpath << ", " << treepath << std::endl;
    int max_events=0; //400000; // limit for development
    try {
        // the categories
#if 0 // gcc objects?
        using ClassifyCal::LOW;
        using ClassifyCal::MED;
        using ClassifyCal::HIGH;
        using ClassifyCal::ALL;
        using ClassifyVertex::THIN;
        using ClassifyVertex::THICK;
        using ClassifyCore::VERTEX;
        using ClassifyCore::TRACK;
#endif
        if( name=="goodcal_low"  || all) ClassifyCal("goodcal_low", ClassifyCal::LOW).run();
        if( name=="goodcal_med"  || all) ClassifyCal("goodcal_med", ClassifyCal::MED).run();
        if( name=="goodcal_high" || all) ClassifyCal("goodcal_high",ClassifyCal::HIGH).run();
//        if( name=="goodcal_all"  || all) ClassifyCal("goodcal_ALL", ALL).run();
        if( name=="goodcal"      || all) ClassifyCal("goodcal").run();

        if( name=="vertex_thin"     || all) ClassifyVertex("vertex_thin",  ClassifyVertex::THIN).run();
        if( name=="vertex_thick"    || all) ClassifyVertex("vertex_thick", ClassifyVertex::THICK).run();

        if( name=="psf_thin_vertex" || all || psf) ClassifyCore("psf_thin_vertex", ClassifyCore::VERTEX, ClassifyVertex::THIN).run(max_events);
        if( name=="psf_thick_vertex"|| all || psf) ClassifyCore("psf_thick_vertex",ClassifyCore::VERTEX, ClassifyVertex::THICK).run();
        if( name=="psf_thin_track"  || all || psf) ClassifyCore("psf_thin_track",  ClassifyCore::TRACK,  ClassifyVertex::THIN).run();
        if( name=="psf_thick_track" || all || psf) ClassifyCore("psf_thick_track", ClassifyCore::TRACK,  ClassifyVertex::THICK).run();

        if( name=="gamma"           || all) ClassifyGamma("gamma").run();
    }
    catch (const std::exception & error)
    {
        std::cerr << "Caught exception \"" << error.what() << "\"" <<std::endl;
        rc= 1;
    }
    return rc;
}
