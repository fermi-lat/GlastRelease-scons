//$Header$
// Original author T. Burnett (w/ help from H. Kelley)
#ifndef ROOTTUPLE_H
#define ROOTTUPLE_H
#include "analysis/Tuple.h"


#include <string>

class TFile;
class TNtuple;
class TTree;

class RootTuple : public Tuple {

public:
    ///
    RootTuple::RootTuple(std::string title, std::string file, std::string treeName);
    ~RootTuple(){};

    //! override this to create from the Gaudi-style n-tuple
    const TupleItem* tupleItem(const std::string& name)const;

    //! return false when no more events
    bool nextEvent();

    //! return false if event idx does not exist
    bool getEvent(int idx);

    int numEvents(){return m_numEvents;}
private:
    TTree * m_tree;

    int m_numEvents;
    int m_event;
};

#endif