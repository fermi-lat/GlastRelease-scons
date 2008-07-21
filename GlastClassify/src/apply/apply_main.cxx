/** @file apply_main.cxx 
@brief Application that applies decision trees to a root tuple

$Header$
*/

#include "RootTuple.h"

#include "../AtwoodTrees.h"
#include "facilities/Util.h"

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

        std::string CTFilePath = "";
        const char* ctree = ::getenv("CTREE_PATH");
        if (ctree==0)
        {
	    //std::string rootPath = ::getenv("GLASTCLASSIFYROOT");
            //facilities::Util::expandEnvVar(&rootPath);

            CTFilePath = facilities::commonUtilities::joinPath(facilities::commonUtiliies::getXmlPath("GlastClassify"), "Pass4_Analysis_Complete.xml");
            std::cout << "Setting file to  " << CTFilePath << std::endl;
        }
        else
        {
            std::cout << "CTREE_PATH: " << ctree << std::endl;
            CTFilePath = std::string(ctree);
        }

        // Expand the environment variable name
        facilities::Util::expandEnvVar(&CTFilePath);

        // create the ct: pass in the tuple.
        AtwoodTrees ctrees(tuple, std::cout, CTFilePath);

        // Are we pruning as well?
        bool keepAllRows = true; // Eventually provide ability to set this true/false...

        const char* pruneRows = ::getenv("PRUNEROWS");
        if (pruneRows)
        {
            std::string pruneEm(pruneRows);

            if (pruneEm == "true") keepAllRows = false;
        }

        // set up the output Root file, branch
        tuple.setOutputFile(output_filename);

        int numInputRows  = 0;
        int numOutputRows = 0;
        
        int numIntervals  = 50;
        int total(tuple.numEvents());
        int fraction = total / numIntervals;

        if (fraction < 2) fraction = 2;

        for( int i = 0; i < numIntervals+2; i++) std::cout << " ";
        std::cout << "]\r[";

        while ( tuple.nextEvent() ) 
        {
            // Execute the tree analysis on this tuple row
            // If good result then also write out the new row
            if (ctrees.execute() || keepAllRows)   // fill in the classification (testing here)
            {
                tuple.fill();
                numOutputRows++;
            }

            if( (numInputRows++ % fraction) == 0 ) std::cout << ".";
        }

        std::cout << "]" << std::endl;
        std::cout << "Read  " << numInputRows << " entries" << std::endl;
        std::cout << "Wrote " << numOutputRows << " entries" << std::endl;

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


