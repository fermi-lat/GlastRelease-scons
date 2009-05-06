//$Header$
// Original author T. Burnett (w/ help from H. Kelley)
#ifndef ROOTTUPLE_H
#define ROOTTUPLE_H
#include "GlastClassify/ITupleInterface.h"


#include <string>

class TFile;
class TChain;
class TTree;

class RootTuple : public GlastClassify::ITupleInterface {

public:
    /// 
    RootTuple::RootTuple( std::string inputfile, std::string treeName);
    ~RootTuple();

    //! acccess to an item (interface to the leaf)
    const GlastClassify::Item* getItem(const std::string& name)const;

    //! create new leaf (float only)
    void addItem(const std::string& name, float& value);
    void addItem(const std::string& name, double& value);
    void addItem(const std::string& name, unsigned long long& value);
    void addItem(const std::string& name, char & value);

    //! return false when no more events
    bool nextEvent();

    //! return false if event idx does not exist
    bool getEvent(int idx);

    int numEvents(){return m_numEvents;}

    //! set the output file name: if set, will copy to it
    void setOutputFile(const std::string& outputfilename);

    //! set current event
    void fill();

private:
    TTree*  tree(){return m_tree;}
    TTree*  m_tree;
    TTree*  m_output_tree;
    TFile*  m_file;
    TChain* m_inChain;
    TFile*  m_output_file;

    int m_numEvents;
    int m_event;
};

#endif
