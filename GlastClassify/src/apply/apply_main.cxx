/** @file apply_main.cxx 
@brief Application that applies decision trees to a root tuple

$Header$
*/

#include "RootTuple.h"

#include "GlastClassify/AtwoodTrees.h"

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>

#include "TFile.h"
#include "TTree.h"

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

        std::stringstream title; 
        title << "gen(" << tuple.numEvents() << ")";


        const char* ctree = ::getenv("CTREE_PATH");
        if (ctree==0) ctree = "D:\\common\\ctree\\GlastClassify\\v1r0\\treeinfo";
        // create the ct: pass in the tuple.
        AtwoodTrees ctrees(tuple, std::cout, ctree!=0? std::string(ctree) : "");

        // set up the output Root file, branch

        TFile out_file(output_filename.c_str(), "recreate");
        TTree* out_tree = tuple.tree()->CloneTree();

        int k(0);
        while ( tuple.nextEvent() ) { 
            ctrees.execute();   // fill in the classification (testing here)
            out_file.cd();
            out_tree->Fill();
            ++k;
        }
        out_file.Write();
        std::cout << "Wrote " << k << " entries" << std::endl;

    }catch(std::exception& e){
        std::cerr << "Caught: " << e.what( ) << std::endl;
        std::cerr << "Type: " << typeid( e ).name( ) << std::endl;
        rc=1;
    }catch(...) {
        std::cerr << "Unknown exception from classification " << std::endl;
        rc=2;
    }
    return rc;
}


