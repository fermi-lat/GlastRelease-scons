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


GlastClassify::GlastClassify(const std::string& info_path, bool mixed)
: m_info(s_treepath+"/"+info_path, s_rootpath)
, m_mixed(mixed)
{
    std::ofstream* logfile = new std::ofstream(m_info.log().c_str());
    if( ! logfile->is_open() ) {
        throw std::runtime_error("log file "+  m_info.log() + " did not open");
    }
    m_nobkgnd = true;
    m_log = logfile;
    std::ifstream titlefile((m_info.filepath()+"/title.txt").c_str());
    titlefile >> m_title;
    log() << "Starting classification of " << m_title << std::endl;  
    std::cout << "Starting classification of " << m_title << std::endl;  
}

void GlastClassify::run( unsigned int max_events, Subset set)
{
    current_time();
    current_time(log());

    // load  the events
    load(max_events, set);

    current_time();
    std::cout << "begin classification" << std::endl;

    // run the classification

    classify();
    current_time();
    current_time(log());
}

int GlastClassify::subdefine(std::vector<std::string>& all_names, const char *Filename)
{
    all_names.push_back(Filename);
    return all_names.size()-1;
}

void GlastClassify::load( unsigned int max_events, Subset /*set*/)
{
    std::vector<std::string> all_names(m_info.vars());
    std::cout << "Defining names" << std::endl;
    define(all_names);
    if(m_mixed){
        // signal and background are mixed in one batch of files
        load(m_info.signalFiles(), all_names);
    }else {
        // separate signal, background
        load(m_info.signalFiles(), all_names, true);
        load(m_info.backgroundFiles(), all_names,  false);
    }
}

void GlastClassify::load(TrainingInfo::StringList input, std::vector<std::string> all_names,  bool isSignal)
{
    std::cout << "loading " << input.size() << " files: \n\t" ;
    std::copy(input.begin(), input.end(), std::ostream_iterator<std::string>(std::cout, "\n\t")); 
    std::cout << std::endl;

    log() << "Processing file(s)\n\t";
    std::copy(input.begin(), input.end(), std::ostream_iterator<std::string>(log(), "\n\t"));

    std::cout << "calling RootTuple" << std::endl;
    RootTuple t(input, "MeritTuple");
    t.selectColumns(all_names, false); // not weighted
    double good=0, bad=0, rejected=0;
    Classifier::Record::setup();
    RootTuple::Iterator rit = t.begin();
    log() << "\tsize = " << t.size() << std::endl;

    unsigned int max_events=0;
    int nvars = all_names.size();
    for( ; rit!=t.end();  ++rit)
    { 
        if (max_events>0 && rit > max_events) break ; //  max
        m_row = &*rit;
        // copy to local
        if( !accept() ) {++rejected; continue; } // apply general cut
        bool signal = m_mixed ? isgood() : isSignal;
        if (signal) ++good; else ++bad;
        // copy the data to the Classifier's table for the classification
        m_data.push_back( Classifier::Record(signal, m_row->begin(), m_row->begin()+nvars));
    }
    log() << "\tgood, bad, rejected records: " << good << ",  " << bad <<", " << rejected << std::endl;
    log() << "Loaded " << (good+bad) << " records"<< std::endl;
}


void GlastClassify::classify()
{
    // create the tree from the data
    Classifier tree(m_data, m_info.vars());
    tree.makeTree();

    // print the node list, and the variables used
    tree.printTree(log());
    tree.printVariables(log());

    // a table of the background for a given efficiency
    BackgroundVsEfficiency plot(tree);
    plot.print(log());
    log() << "Figure of merit sigma: " <<  plot.sigma() << std::endl;
    std::ofstream plotfile((m_info.filepath()+"/efficiencyplot.txt").c_str());
    plot.print(plotfile);
    

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
    if( x== b ) throw std::runtime_error(std::string("could not find variable ")+name);
    return x  - a;
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



