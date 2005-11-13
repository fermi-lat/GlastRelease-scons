/** @file apply_main.cxx 
@brief Application that applies decision trees to a root tuple

$Header$
*/

#include "RootTuple.h"

#include "GlastClassify/AtwoodLikeTrees.h"

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>

/** 
@page applications apply
format: apply [input_file] [output_file]

@param input_file if not present look at env var MERIT_INPUT_FILE
@param output_file if not present, and MERIT_OUTPUT_FILE is not defined, just append "_new" to the file name.

Copies the root tree MeritTuple from input_file to output_file, but recalculates the CT variables,
Expects the env var CTREE_PATH to point to a folder containing the trees.


*/

int main(int argc, char* argv[])
{
    using namespace GlastClassify;

    int rc = 0;
    try {

        std::string  input_filename(""), output_filename(""), tree_name("MeritTuple");
        int n=0;
        if( argc>++n ) input_filename = argv[n];		// required
        if( argc>++n ) output_filename = argv[n];		// required

        if( input_filename=="" ) {
            const char * env = ::getenv("MERIT_INPUT_FILE");
            if( env ) input_filename=env;
            else {
                throw std::invalid_argument( "No input file specified");
            }
        }
        if( output_filename=="" ) {
            const char * env = ::getenv("MERIT_OUTPUT_FILE");
            if( env ) output_filename=env;
            else {
                // make up output file from input
                int find = input_filename.find(".root");
                output_filename = input_filename.substr(0,find)+"_new.root";
            }
        }

        std::cerr << "Converting CT variables from: \"" << input_filename << "\" to\n\t\"" 
            << output_filename << "\"" << std::endl;

        RootTuple tuple(input_filename, tree_name);

        const char* ctree = ::getenv("CTREE_PATH");
        if (ctree==0){
            ctree = "D:\\common\\ctree\\GlastClassify\\v2r1\\treeinfo-v4";
            std::cout << "Setting CTREE_PATH to  " << ctree << std::endl;
        }else{
            std::cout << "CTREE_PATH: " << ctree << std::endl;
            }

        // create the ct: pass in the tuple.
        AtwoodLikeTrees ctrees(tuple, std::cout, std::string(ctree) );

        // set up the output Root file, branch

        tuple.setOutputFile(output_filename);
        int k(0), total(tuple.numEvents()), fraction(0);
        for( int i = 0; i<51; ++i)std::cout << " ";
        std::cout << "]\r[";
        while ( tuple.nextEvent() ) { 
            ctrees.execute();   // fill in the classification (testing here)
            tuple.fill();
            ++k;
            if(  (50*k)% total == 0 ) {
                std::cout << ".";
            }
        }
        std::cout << "]\nWrote " << k << " entries" << std::endl;

    }catch(const std::exception& e){
        std::cerr << "Caught: " << e.what( ) << std::endl;
        std::cerr << "Type: " << typeid( e ).name( ) << std::endl;
        rc=1;
    }catch(...) {
        std::cerr << "Unknown exception from classification " << std::endl;
        rc=2;
    }    
    return rc;
}


