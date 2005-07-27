/**@file main.cxx
@brief main program for application to create GLAST classification trees
$Heading$
*/

#include "ClassifyCal.h"
#include "ClassifyVertex.h"
#include "ClassifyCore.h"
#include "ClassifyGamma.h"

int main(int argc , char * argv[])
{
    int rc=0;
    std::string // defaults
        rootpath("d:\\common\\DC2\\root_files"),
        treepath("d:\\common\\ctree\\data"),
        name("goodcal");
    if( argc>1) rootpath=argv[1];
    if( argc>2) treepath=argv[2];
    if( argc>3) name = argv[3];
    bool all =  name=="all"; 

    GlastClassify::setPaths(rootpath, treepath);

    std::cout << "classifier invoked with root, tree paths: "<< rootpath << ", " << treepath << std::endl;
    try {
        if( name=="goodcal"         || all ) ClassifyCal("goodcal").run();

        if( name=="vertex_thin"     || all)  ClassifyVertex("vertex_thin",  ClassifyVertex::layer::THIN).run();
        if( name=="vertex_thick"    || all)  ClassifyVertex("vertex_thick", ClassifyVertex::layer::THICK).run();

        if( name=="psf_thin_vertex" || all) ClassifyCore("psf_thin_vertex", ClassifyCore::etype::VERTEX, ClassifyVertex::layer::THIN).run();
        if( name=="psf_thick_vertex"|| all) ClassifyCore("psf_thick_vertex",ClassifyCore::etype::VERTEX, ClassifyVertex::layer::THICK).run();
        if( name=="psf_thin_track"  || all) ClassifyCore("psf_thin_track",  ClassifyCore::etype::TRACK, ClassifyVertex::layer::THIN).run();
        if( name=="psf_thick_track" || all) ClassifyCore("psf_thick_track", ClassifyCore::etype::TRACK, ClassifyVertex::layer::THICK).run();
        
        if( name=="gamma"           || all) ClassifyGamma("gamma").run();
    }
    catch (const std::exception & error)
    {
        std::cerr << "Caught exception \"" << error.what() << "\"" <<std::endl;
        rc= 1;
    }
    return rc;
}
