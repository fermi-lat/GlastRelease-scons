/**@file GlastClassify.cxx
@brief 

$Header$
*/

#include "GlastClassify.h"

#include "classifier/DecisionTree.h"
#include "classifier/BackgroundVsEfficiency.h"
#include "classifier/TrainingInfo.h"
#include "classifier/RootTuple.h"
#include "classifier/Trainer.h"
#include "classifier/AdaBoost.h"

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
int GlastClassify::s_boost=0;
int GlastClassify::s_events=0;
int GlastClassify::s_normalize=100;

GlastClassify::Subset GlastClassify::s_train_sample=ODD;
GlastClassify::Subset GlastClassify::s_test_sample=EVEN;

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
, m_total_good(0), m_total_bad(0)
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
        <<"Starting "<<(s_train? "classification and testing ":"testing only of " )<< m_title << std::endl;  
    std::cout << "Starting "<<(s_train? "classification and testing ":"testing of " )<< m_title << std::endl;  
}

void GlastClassify::run()
{
    current_time();
    current_time(log());

    // run the train or test (or both)?

    if( s_train) {
        load(s_events, s_train_sample);
        classify();
    }
    log() << "========================== testing ============================\n";
    std::cout << "start testing ...\n";
    load( s_events, s_test_sample);
    test();

    current_time();
    current_time(log());
}

int GlastClassify::add_index(const std::string& name)
{
    m_all_names.push_back(name);
    return m_all_names.size()-1;
}

void GlastClassify::load( unsigned int max_events , Subset set)
{
    log() << "Loading ";
    switch (set){
        case ALL:  log() << "all"; break;
        case EVEN: log() << "even"; break;
        case ODD : log() << "odd"; break;
            //       case RANDOM : log() << "random"; break;
    }
    log() << " events." << std::endl;

    m_data.clear();
    m_total_good= m_total_bad=0;

    if(m_mixed){
        // signal and background are mixed in one batch of files
        load(m_info.signalFiles(), max_events, set);
    }else {
        // separate signal, background
        load(m_info.signalFiles(),    max_events, set, true);
        load(m_info.backgroundFiles(),max_events, set, false);
    }
    log() << "Event totals: good "<< m_total_good << ", bad " << m_total_bad << std::endl;

    if( s_normalize>0) {
        log() << "Normalizing signal, background to be equal to "<< s_normalize << std::endl;
        m_data.normalize(s_normalize, s_normalize);

    }
}
void GlastClassify::load(TrainingInfo::StringList input, 
                         unsigned int max_events,  Subset set, bool isSignal)
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
    log() << "\tsize = " << t.size() << std::endl;
    int nvars = m_all_names.size();

    RootTuple::Iterator rit = t.begin();
    if( set==EVEN) ++rit; // skip the first if EVEN
    //if( set==RANDOM && m_rand->shoot()>0.5) ++rit; // skip first

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
        if( accept() ) {
           // copy to local
            bool signal = m_mixed ? isgood() : isSignal;
            if (signal) ++good; else ++bad;
            // copy the data to the Classifier's table for the classification
            m_data.push_back( Classifier::Record(signal, m_row->begin(), m_row->begin()+nvars));
        } else {
           ++rejected;
         }
        // if doing alternate or random, skip the next record here.
        if( set !=ALL ){  //   || set==RANDOM && m_rand->shoot()>0.5 
            ++rit; if( rit >= t.end() )break;
        }
    }
    log() << "\tgood, bad, rejected records: " << good << ",  " << bad <<", " << rejected << std::endl;
    if( nan>0) log() << "\tWARNING: found "<< nan << " events with non-finite values " << std::endl;
    log() << "Loaded " << (good+bad) << " records"<< std::endl;
    m_total_good += good;
    m_total_bad += bad;

}


void GlastClassify::classify()
{

    // create the tree from the data
    Classifier ctree(m_data, m_info.vars());
    ctree.makeTree();

    // summary stuff at the top of the file
    log() << "Number of nodes in the tree: " << Classifier::Node::s_nodes << std::endl;

    // insert list of used variables in the log, and to a file
    ctree.printVariables(log());
    log() << "======================================\n";
    std::string usedvarfilename(m_info.filepath()+"/used_variables.txt");
    std::ofstream usedvarsfile(usedvarfilename.c_str());
    ctree.printVariables(usedvarsfile);
    log()<<" writing to file " << usedvarfilename << std::endl;


    BackgroundVsEfficiency plot(ctree);
    log() << "\nFigure of merit sigma: " <<  plot.sigma() << std::endl;

    std::string plotfilename(m_info.filepath()+"/train_efficiency.txt");
    log()<<" writing to file " << plotfilename << std::endl;

    std::ofstream plotfile(plotfilename.c_str());
    plot.print(plotfile);

#ifdef VERBOSE // generates a lot of output, need a special option
    // print the node list, and the variables used
    tree.printTree(log());
#endif

    if( s_boost==0) {

        // single tree, no boosting
        m_dtree = ctree.createTree(m_info.title());

    }else{
        boost(ctree);
        BackgroundVsEfficiency plot(ctree);

        // redo, and overwrite, the plot
        log() << "\nFigure of merit after boosting: " <<  plot.sigma() << std::endl;
        std::string plotfilename(m_info.filepath()+"/train_efficiency.txt");
        log()<<" writing to file " << plotfilename << std::endl;

        std::ofstream plotfile(plotfilename.c_str());
        plot.print(plotfile);

    }
    std::string dtfilename(m_info.filepath()+"/dtree.txt");
    std::ofstream dtfile(dtfilename.c_str());
    log()<<" writing tree to file " << dtfilename << std::endl;
    m_dtree->print(dtfile);
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

void GlastClassify::boost(Classifier & classify)
{
    // set up the booster, and create the first boosted tree
    //
    AdaBoost booster(m_data, 0.5); //TODO adaBeta);
    double boostwt = booster(classify); //weight of initial tree, boost training sample
    m_dtree = classify.createTree(m_info.title(), boostwt); 

    for (int itree = 1; itree < s_boost+1; ++itree) {
        if( itree%5==1 ){ // checkpoint every 5, starting with single tree
            //TODO 
#if 0
            saveTree();
            double resolution = evaluate(testing_sample) ;
            std::cout << "Test resolution: " << resolution << std::endl;
            log() << "Test resolution: " << resolution << std::endl;
#endif
        }
        log() << "Making boosted tree #" 
            << itree << ", current weight " << boostwt  
            << ", nodes: "<<Classifier::Node::s_nodes 
            << std::endl;
        std::cout << "Making boosted tree #" 
            << itree << ", current weight " << boostwt 
            << ", nodes: " << Classifier::Node::s_nodes 
            << std::endl;
        if ( boostwt < 1.001) {
            log() << "quitting boost: weight too small to continue" << std::endl;
            break;
        }
        // create new Classifier object with reweighted data from booster
        Classifier classify(booster.data());
        classify.makeTree(); // and make a new tree using it

        // get the weight from the booster, which also reweights data for next cycle
        boostwt = booster(classify);

        // generate the decision tree from the current classification, append it 
        std::auto_ptr<DecisionTree> boostedtree( classify.createTree(m_info.title(),boostwt));
        m_dtree->addTree(&*boostedtree); //(how to pass a pointer to an auto_ptr object)

    }
}

