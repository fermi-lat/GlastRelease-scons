/**@file main.cxx
@brief main program for application to create GLAST classification trees

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
        treepath("d:\\common\\ctree\\data");
    if( argc>1) rootpath=argv[1];
    if( argc>2) treepath=argv[2];

    GlastClassify::s_rootpath = rootpath+"/";
    GlastClassify::s_treepath = treepath+"/";
    try {
#if 0
        ClassifyCal(cpath+"goodcal").run();
        ClassifyVertex("../data/vertex_thin", ClassifyVertex::layer::THIN).run();
        ClassifyVertex("../data/vertex_thick", ClassifyVertex::layer::THICK).run();
        ClassifyCore("../data/psf_thin_vertex", ClassifyCore::etype::VERTEX, ClassifyVertex::layer::THIN).run();
        ClassifyCore("../data/psf_thick_vertex", ClassifyCore::etype::VERTEX, ClassifyVertex::layer::THICK).run();
        ClassifyCore("../data/psf_thin_track", ClassifyCore::etype::TRACK, ClassifyVertex::layer::THIN).run();
        ClassifyCore("../data/psf_thick_track", ClassifyCore::etype::TRACK, ClassifyVertex::layer::THICK).run();
#else		
        ClassifyGamma("gamma").run();
#endif
    }
	catch (const std::exception & error)
	{
        std::cerr << "Caught exception \"" << error.what() << "\"" <<std::endl;
        rc= 1;
    }
    return rc;
}
