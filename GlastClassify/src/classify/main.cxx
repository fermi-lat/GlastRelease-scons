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
        rootpath("F:\\glast\\data"),
        treepath("..\\data"),
        name("all");
    bool train(true);
    if( argc==1) {
        // no args: check env vars
        const char* env = ::getenv("ROOTPATH");
        if (env!=0) rootpath = std::string(env);
        env = ::getenv("TREEPATH");
        if( env!=0) treepath = std::string(env);
        env = ::getenv("CLASSIFY_TYPE");
        if( env!=0) name = std::string(env);
    }else {
        if( argc>1 && std::string(argv[1]) !="-") rootpath=argv[1];
        if( argc>2 && std::string(argv[2]) !="-") treepath=argv[2];
        if( argc>3) name   = argv[3];
        if( argc>4) GlastClassify::s_train = false; //any char triggers test:assume tree exists
    }
    bool all(name=="all"); 
    bool psf(name=="psf");
    bool energy(name=="energy");
    bool vertex(name=="vertex");

    GlastClassify::setPaths(rootpath, treepath);

    std::cout << "classifier invoked with root, tree paths: "<< rootpath << ", " << treepath << std::endl;
    std::cout << " \tIn "<< (GlastClassify::s_train? " training" : "testing") << " mode" << std::endl;
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
        if( name=="goodcal_low"  || all || energy) ClassifyCal("goodcal_low", ClassifyCal::LOW).run();
        if( name=="goodcal_med"  || all || energy) ClassifyCal("goodcal_med", ClassifyCal::MED).run();
        if( name=="goodcal_high" || all || energy) ClassifyCal("goodcal_high",ClassifyCal::HIGH).run();

        if( name=="vertex_thin"     || all|| vertex) ClassifyVertex("vertex_thin",  ClassifyVertex::THIN).run();
        if( name=="vertex_thick"    || all|| vertex) ClassifyVertex("vertex_thick", ClassifyVertex::THICK).run();

        if( name=="psf_thin_vertex" || all || psf) ClassifyCore("psf_thin_vertex", ClassifyCore::VERTEX, ClassifyVertex::THIN).run(max_events);
        if( name=="psf_thick_vertex"|| all || psf) ClassifyCore("psf_thick_vertex",ClassifyCore::VERTEX, ClassifyVertex::THICK).run();
        if( name=="psf_thin_track"  || all || psf) ClassifyCore("psf_thin_track",  ClassifyCore::TRACK,  ClassifyVertex::THIN).run();
        if( name=="psf_thick_track" || all || psf) ClassifyCore("psf_thick_track", ClassifyCore::TRACK,  ClassifyVertex::THICK).run();

        if( name.substr(0,6)=="gamma/" || all) ClassifyGamma(name).run();
        if( name=="gamma" || all ) {
            std::string gamma_subset_name[]={"highcal", "medcal", "thin", "thick"};
            for( int i=0; i<4; ++i){
                ClassifyGamma("gamma/vertex/"+gamma_subset_name[i]).run();
                ClassifyGamma("gamma/track/"+gamma_subset_name[i]).run();
            }
        }
    }
    catch (const std::exception & error)
    {
        std::cerr << "Caught exception \"" << error.what() << "\"" <<std::endl;
        rc= 1;
    }
    return rc;
}
