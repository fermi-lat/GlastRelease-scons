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


*/
int main(int argc , char * argv[])
{
    int rc=0;
    std::string // defaults for development
        rootpath("d:\\common\\DC2\\root_files"),
        treepath("d:\\common\\ctree\\data"),
        name("all");
    if( argc>1) rootpath=argv[1];
    if( argc>2) treepath=argv[2];
    if( argc>3) name = argv[3];
    bool all =  name=="all"; 

    GlastClassify::setPaths(rootpath, treepath);

    std::cout << "classifier invoked with root, tree paths: "<< rootpath << ", " << treepath << std::endl;
    try {
        // the categories
        using ClassifyCal::LOW;
        using ClassifyCal::MED;
        using ClassifyCal::HIGH;
        using ClassifyCal::ALL;
        using ClassifyVertex::THIN;
        using ClassifyVertex::THICK;
        using ClassifyCore::VERTEX;
        using ClassifyCore::TRACK;

        if( name=="goodcal_low"  || all) ClassifyCal("goodcal_low", LOW).run();
        if( name=="goodcal_med"  || all) ClassifyCal("goodcal_med", MED).run();
        if( name=="goodcal_high" || all) ClassifyCal("goodcal_high",HIGH).run();
        if( name=="goodcal_all"  || all) ClassifyCal("goodcal_ALL", ALL).run();
        if( name=="goodcal"      || all) ClassifyCal("goodcal").run();

        if( name=="vertex_thin"     || all) ClassifyVertex("vertex_thin",  THIN).run();
        if( name=="vertex_thick"    || all) ClassifyVertex("vertex_thick", THICK).run();

        if( name=="psf_thin_vertex" || all) ClassifyCore("psf_thin_vertex", VERTEX, THIN).run();
        if( name=="psf_thick_vertex"|| all) ClassifyCore("psf_thick_vertex",VERTEX, THICK).run();
        if( name=="psf_thin_track"  || all) ClassifyCore("psf_thin_track",  TRACK,  THIN).run();
        if( name=="psf_thick_track" || all) ClassifyCore("psf_thick_track", TRACK,  THICK).run();
        if( name=="gamma"           || all) ClassifyGamma("gamma").run();
    }
    catch (const std::exception & error)
    {
        std::cerr << "Caught exception \"" << error.what() << "\"" <<std::endl;
        rc= 1;
    }
    return rc;
}
