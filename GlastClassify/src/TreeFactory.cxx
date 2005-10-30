/**@file TreeFactory.cxx

@brief implementation of class TreeFactory

$Header$
*/

#include "GlastClassify/TreeFactory.h"
#include "classifier/DecisionTree.h"
#include "classifier/TrainingInfo.h"
#include "Classifier/Filter.h"

#include <fstream>
#include <cassert>
#include <stdexcept>

using namespace GlastClassify;

/** @class GleamValues
@brief local definition of class which handles pointers to values
*/
class TreeFactory::GleamValues : public DecisionTree::Values {
public:
    GleamValues(const TrainingInfo::StringList & names,ILookupData& lookup)
    {
        for( TrainingInfo::StringList::const_iterator it = names.begin();
            it != names.end(); ++it)
        {
            std::pair<bool, const void*> entry = lookup(*it);
            if( entry.second == 0 ){
                throw std::invalid_argument("TreeFactory::GleamValues: did not find variable "+*it);
            }

            m_pval.push_back( entry);
        }

    }
    /// @brief callback from tree evaluation

    double operator[](int index)const{
        std::pair<bool, const void*> entry = m_pval[index];
        // now dereference either as a float or a double
        double v = entry.first? *(const float*)entry.second : *(const double*)entry.second;
        return v;
    }
private:
    std::vector<std::pair<bool, const void*> >m_pval;
};


const TreeFactory::Tree& TreeFactory::operator()(const std::string& name)
{
    m_trees.push_back(new Tree(m_path+"/"+name, m_lookup));
    return *m_trees.back();
}

TreeFactory::Tree::Tree( const std::string& path, ILookupData& lookup)
{
    TrainingInfo info(path);
    TrainingInfo::StringList vars(info.vars()); // local copy to extend perhaps

    // make sure can open the classification tree file
    std::string dtfile(path+"/dtree.txt");
    std::ifstream dtstream(dtfile.c_str());
    if(! dtstream.is_open()){
        throw std::runtime_error("failed to open file "+dtfile);
    }
    DecisionTree* theTree = new DecisionTree(dtstream);

    // see if there is a filter

    std::string filterfile(path+"/filter.txt");
    std::ifstream filterstream(filterfile.c_str());
    if( filterstream.is_open()){
        filterstream.close();

        // yes, create tree object with it-- note that vars may be extended   
        DecisionTree* dt = new DecisionTree(info.title());
        Filter filter(vars,*dt); 
        filter.addCutsFrom(filterfile);
        filter.close();

        // and add the rest here
        dt->addTree(theTree);
        delete theTree; // done with this one
        m_dt=dt;        // copy pointer: now it is const.

    }else{

        // no filter: this is it.
        m_dt=theTree;
    }
    m_vals = new GleamValues(vars, lookup);
}

double TreeFactory::Tree::operator()()const
{
    return (*m_dt)(*m_vals);
}
std::string TreeFactory::Tree::title()const
{
    return m_dt->title();
}
TreeFactory::Tree::~Tree()
{
    delete m_dt;
    delete m_vals;
}
TreeFactory::~TreeFactory()
{
    ///@TODO: delete trees
}
