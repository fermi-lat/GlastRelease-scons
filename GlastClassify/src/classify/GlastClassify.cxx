/**@file GlastClassify.cxx
@brief 

$Header$
*/

#include "GlastClassify.h"

#include "classifier/DecisionTree.h"
#include "classifier/BackgroundVsEfficiency.h"
#include "classifier/TrainingInfo.h"
#include "classifier/RootTuple.h"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <time.h>
#include <stdio.h>

    /// path to the root files
std::string GlastClassify::s_rootpath;
    /// path to tree data, input and output
std::string GlastClassify::s_treepath;

bool GlastClassify::s_train(true); /// set false for test mode.

namespace {
    GlastClassify* s_current;
}


GlastClassify::Entry::Entry( const std::string& name)
: m_cl(s_current) // avoid warning with kluge
, m_index(m_cl->find_index(name))
{
}




GlastClassify::GlastClassify(const std::string& info_path, bool mixed)
: m_info(s_treepath+"/"+info_path, s_rootpath)
, m_all_names(m_info.vars()) 
, m_mixed(mixed)
{
    s_current= this; // for communication of pointer during setup

    std::ofstream* logfile = new std::ofstream(m_info.log().c_str(), std::ofstream::app);
    if( ! logfile->is_open() ) {
        throw std::runtime_error("log file "+  m_info.log() + " did not open");
    }
    m_nobkgnd = true;
    m_log = logfile;
    std::ifstream titlefile((m_info.filepath()+"/title.txt").c_str());
    titlefile >> m_title;
    log() << "=======================================================================\n"
        <<"Starting "<<(s_train? "classification ":"testing of " )<< m_title << std::endl;  
    std::cout << "Starting "<<(s_train? "classification ":"testng of " )<< m_title << std::endl;  
}

void GlastClassify::run( unsigned int max_events, Subset set)
{
    current_time();
    current_time(log());

    // load  the events
    load(max_events, set);

    current_time();

    // run the train or test

    if( s_train) {
        classify();
    }else{
        test();
    }
    current_time();
    current_time(log());
}

int GlastClassify::add_index(const std::string& name)
{
    m_all_names.push_back(name);
    return m_all_names.size()-1;
}

void GlastClassify::load( unsigned int max_events , Subset /*set*/)
{
    if(m_mixed){
        // signal and background are mixed in one batch of files
        load(m_info.signalFiles(), max_events);
    }else {
        // separate signal, background
        load(m_info.signalFiles(), max_events, true);
        load(m_info.backgroundFiles(),max_events, false);
    }
}

void GlastClassify::load(TrainingInfo::StringList input, unsigned int max_events,  bool isSignal)
{
    std::cout << "loading " << input.size() << " files: \n\t" ;
    std::copy(input.begin(), input.end(), std::ostream_iterator<std::string>(std::cout, "\n\t")); 
    std::cout << std::endl;

    log() << "Processing file(s)\n\t";
    std::copy(input.begin(), input.end(), std::ostream_iterator<std::string>(log(), "\n\t"));

    RootTuple t(input, "MeritTuple");
    t.selectColumns(m_all_names, false); // not weighted
    double good=0, bad=0, rejected=0, nan=0;
    Classifier::Record::setup();
    RootTuple::Iterator rit = t.begin();
    log() << "\tsize = " << t.size() << std::endl;

    int nvars = m_all_names.size();
    for( ; rit!=t.end();  ++rit)
    { 
        if (max_events>0 && rit > max_events) break ; //  max
        try {
            m_row = &*rit;
        } catch( std::runtime_error err){
            ++nan;
            ++rejected;
            continue;
        }
        // copy to local
        if( !accept() ) {++rejected; continue; } // apply general cut
        bool signal = m_mixed ? isgood() : isSignal;
        if (signal) ++good; else ++bad;
        // copy the data to the Classifier's table for the classification
        m_data.push_back( Classifier::Record(signal, m_row->begin(), m_row->begin()+nvars));
    }
    log() << "\tgood, bad, rejected records: " << good << ",  " << bad <<", " << rejected << std::endl;
    if( nan>0) log() << "\tWARNING: found "<< nan << " events with non-finite values " << std::endl;
    log() << "Loaded " << (good+bad) << " records"<< std::endl;
}


void GlastClassify::classify()
{

    // create the tree from the data
    Classifier tree(m_data, m_info.vars());
    tree.makeTree();

    // summary stuff at the top of the file
    tree.printVariables(log());
    log() << "======================================\n";
    BackgroundVsEfficiency plot(tree);
    log() << "\nFigure of merit sigma: " <<  plot.sigma() << std::endl;

    std::string plotfilename(m_info.filepath()+"/train_efficiency.txt");
    log()<<" writing to file " << plotfilename << std::endl;

    std::ofstream plotfile(plotfilename.c_str());
    plot.print(plotfile);
   
#if 0 // generates a lot of output, need a special option
    // print the node list, and the variables used
    tree.printTree(log());
#endif

    // make simple tree and print it to a local file
    DecisionTree& dtree = *tree.createTree(m_info.title());
    std::ofstream dtfile((m_info.filepath()+"/dtree.txt").c_str());
    dtree.print(dtfile);
}

int GlastClassify::find_index(const std::string& name)
{
    std::vector<std::string>::const_iterator 
        a = m_info.vars().begin(),
        b = m_info.vars().end(),

        x = std::find(a, b, name);
    if( x!= b ){
        return x  - a;
    } else {
        return add_index(name);
    }
}

void GlastClassify::current_time(std::ostream& out)
{   
    static bool first=true;
    static time_t start;
    if(first){ first=false; ::time(&start);}
    time_t aclock;
    ::time( &aclock );   
    char tbuf[25]; ::strncpy(tbuf, asctime( localtime( &aclock ) ),24);
    tbuf[24]=0;
    out<<  "Current time: " << tbuf
        << " ( "<< ::difftime( aclock, start) <<" s elapsed)" << std::endl;
}

void GlastClassify::test()
{
    // create the tree from by reading the dtree file
    std::ifstream treefile((m_info.filepath()+"/dtree.txt").c_str());
    DecisionTree& dtree = *new DecisionTree(treefile);
    log() << "======================================\n";
    BackgroundVsEfficiency plot(dtree, m_data);
    log() << "\nFigure of merit sigma: " <<  plot.sigma() << std::endl;

    std::string plotfilename(m_info.filepath()+"/test_efficiency.txt");
    log()<<" writing to file " << plotfilename << std::endl;

    std::ofstream plotfile(plotfilename.c_str());

    // a table of the background for a given efficiency
    plot.print(plotfile);

}


