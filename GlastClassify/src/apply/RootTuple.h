//$Header$
// Original author T. Burnett (w/ help from H. Kelley)
#ifndef ROOTTUPLE_H
#define ROOTTUPLE_H
#include "GlastClassify/ITupleInterface.h"


#include <string>

class TFile;
class TNtuple;
class TTree;

class RootTuple : public GlastClassify::ITupleInterface {

public:
    ///
    RootTuple::RootTuple(std::string title, std::string file, std::string treeName);
    ~RootTuple(){};

    //! acccess to an item (interface to the leaf)
    const GlastClassify::Item* getItem(const std::string& name)const;

    //! create new leaf
    void addItem(const std::string& name, float& value);

    //! return false when no more events
    bool nextEvent();

    //! return false if event idx does not exist
    bool getEvent(int idx);

    int numEvents(){return m_numEvents;}

    TTree * tree(){return m_tree;}
private:
    TTree * m_tree;
    TFile * m_file;

    int m_numEvents;
    int m_event;
};

#endif
