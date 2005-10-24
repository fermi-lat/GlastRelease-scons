/** @file ParseOptions.h


*/

#include "GlastClassify.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

/** @page applications classify
Run one or more classifications from MeritTuple root files

classify datapath [treepath]  [--case case] [--train | --test] [--normalize x] [--boost nboost]
              [--leaf-purity-limit purity]

@param datapath path to the root data
@param treepath path to tree parameters, place to put results [default ../infopath].
@param case One of: all, vertex, psf, gamma ... [default all]


*/

class ParseOptions {
public:
    
    int n_boosts;
    int events;
    int normalize;
    std::string datapath, infopath, casename;
    std::vector<std::string> file_paths;
    bool train;
    bool valid;
    double leaf_purity ;
    double improv_min;
    GlastClassify::Subset train_sample;
    GlastClassify::Subset test_sample;
#if 0 // not implemented yet
    int minleafsize;
    std::string merit;
    double ada_boost ;
    std::string testfiles;
#endif


    ParseOptions(int argc, char* argv[])
        : n_boosts( 0)
        , train(true)
        , events(0)
        , normalize(100)
        , valid(false)
        , casename("all")
        , infopath("../treeinfo")
        , leaf_purity(0.5)
        , improv_min(0.0)
        , train_sample(GlastClassify::EVEN)
        , test_sample(GlastClassify::ODD)
#if 0 // not implemented yet
        , minleafsize(100)
        , merit("gini")
        , ada_boost( 0.5)
        , testfiles ("files.txt")
#endif 
 
    {
        for (int iarg = 1; iarg < argc; iarg++) {
            std::string arg = argv[iarg]; 
            if        (arg == "--boost"            ) { atof(argv[++iarg], n_boosts);
            } else if (arg == "--case"             ) { casename = argv[++iarg];
            } else if (arg == "--normalize"        ) { atof(argv[++iarg], normalize);
            } else if (arg == "--events"           ) { atof(argv[++iarg], events);
            } else if (arg == "--test"             ) { train = false;
            } else if (arg == "--leaf-purity-limit") { atof(argv[++iarg], leaf_purity);
            } else if (arg == "--improv-min"       ) { atof(argv[++iarg], improv_min);
            } else if (arg == "--train-sample"     ) { train_sample = parse_subset (argv[++iarg]);
            } else if (arg == "--test-sample"      ) { test_sample = parse_subset (argv[++iarg]);
#if 0 // not yet implemented
            } else if (arg == "--merit"            ) { merit = argv[++iarg];
            } else if (arg == "--min-leaf-size"    ) { atof(argv[++iarg], minleafsize);
            } else if (arg == "--ada-boost"        ) { atof(argv[++iarg], ada_boost); 
            } else if (arg == "--train-sample"     ) { train_sample = parse_subset (argv[++iarg]);
            } else if (arg == "--test-files"       ) { testfiles = argv[++iarg];
#endif
            } else if (arg == "--help" || arg == "-h" 
                                     || arg == "-?") { help(); return;
            } else if (arg[0]=='-'                 ) { std::cout << "qualifier \"" << arg <<"\" not recognized!" << std::endl;
                                                      help(); return;
            } else {                                 file_paths.push_back(arg);
            }
        }
        ///
        /// Do argument consistency checks
        ///

        if( file_paths.size() >0) datapath = file_paths[0];
        if( file_paths.size() >1) infopath = file_paths[1];
#if 0
        if (merit != "gini" && merit != "entropy") {
            std::cout << "ERROR: You must specify the figure of merit as 'gini' or 'entropy'! '" + merit + "' is illegal" 
                << std::endl;
            help(); return;
        }
#endif
        valid=true;

    }
    template <typename X>
    void atof(const char * arg, X& output)
    {
        std::istringstream input(arg); 
        input >> output;
    }

    void help()
    {
        std::cout <<
            "Usage: classify [options] datapath outputpath\n"
            "\n"
            "   datapath -- path to data files\n"
            "   outputpath -- path to folder containing cases.txt\n"
            " \n"
            " Options: \n"
            "  --case name            Case to try, defaults to all.\n"
            "  --boost n              Number of boosts to try, if zero or left off no\n"
            "                          boosting is done\n"
            "  --events  n            Maximum number of records to examine\n"
            "  --normalize x          Normalization factor to apply to both singal, noise. \n"
            "  --leaf-purity-limit n  The leaf purity limit used during boosting (0.5)\n"
            "  --improv-min x         Minimum improvement for splitting a node (0), fraction of total gini, or sum of weights\n"

#if 0 // not yet implemented
            "  --min-leaf-size n      Minimum # of events that can be in a leaf (not weight)\n"
            "                          Defaults to 100\n"
            "  --merit [gini|entropy] What to use to find the cut value.\n"
            "  --ada-boost n          Set the beta boost parameter (0.5 default)\n"
            "  --test-files <fname>   A file containing a list of files to test training against\n"
            "                         Defaults to files.txt\n"
#endif
            ;
    }


    GlastClassify::Subset parse_subset (const std::string &s)
    {
        if (       s == "ALL") {		return GlastClassify::ALL;
        } else if (s == "ODD") {	return GlastClassify::ODD;
        } else if (s == "EVEN") {	return GlastClassify::EVEN;
#if 0
        } else if (s == "FIRSTHALF") {	return GlastClassify::FIRSTHALF;
        } else if (s == "SECONDHALF") {	return GlastClassify::SECONDHALF;
        } else if (s == "RANDOM") {	return GlastClassify::RANDOM;
#endif
        } else {
            std::cout << "Illegal subset '" << s << "' -- must be ALL, ODD, EVEN, FIRSTHALF, SECONDHALF, or RANDOM" << std::endl;
            throw std::runtime_error ("Bad argument for subset!");
        }
    }


};


