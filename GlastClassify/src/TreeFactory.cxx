/**@file TreeFactory.cxx

@brief implementation of class TreeFactory

$Header$
*/

#include "GlastClassify/TreeFactory.h"
#include "classifier/DecisionTree.h"
#include "classifier/TrainingInfo.h"
#include "classifier/Filter.h"

#include <fstream>
#include <cassert>
#include <stdexcept>

using namespace GlastClassify;

/** @class GleamValues
@brief local definition of class which handles pointers to values
*/
class TreeFactory::GleamValues : public DecisionTree::Values {
public:
    GleamValues(const TrainingInfo::StringList & names, ITupleInterface& tuple)
    {
        for( TrainingInfo::StringList::const_iterator it = names.begin();
            it != names.end(); ++it)
        {
            const Item * item = tuple.getItem(*it);
            if( item==0){
                throw std::invalid_argument("TreeFactory::GleamValues: did not find variable "+*it);
            }

            m_pval.push_back( item);
        }

    }
    /// @brief callback from tree evaluation

    double operator[](int index)const{

        return *(m_pval[index]);
    }
private:
    std::vector<const Item *  >m_pval;
};


double TreeFactory::evaluate(int i)const {return (*m_trees[i])();}

#if 0
const ITreeFactory::ITree& TreeFactory::operator()(const std::string& name)
#else
const TreeFactory::Tree& TreeFactory::operator()(const std::string& name)

#endif
{
    m_trees.push_back(new Tree(m_path+"/"+name, m_tuple));
    return *m_trees.back();
}

TreeFactory::Tree::Tree( const std::string& path, ITupleInterface& tuple)
: m_filter_tree()
{
    TrainingInfo info(path);
    TrainingInfo::StringList vars(info.vars()); // local copy to extend perhaps

    // see if we have a local  classification tree file
    DecisionTree* localTree(0);
    std::string dtfile(path+"/dtree.txt");
    std::ifstream dtstream(dtfile.c_str());
    if(! dtstream.is_open()){

        //throw std::runtime_error("failed to open file "+dtfile);
    }else{
        localTree = new DecisionTree(dtstream);
    }

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
        
        if( localTree!=0){

            // there is a local tree: append it to the filter
            dt->addTree(localTree);
            delete localTree; // done with this one
            m_filter_tree = dt;
        }else {
            // no local tree, so only filter.
           m_filter_tree=dt;        
        }

    }else{

        // no filter: this is it.
        m_filter_tree = localTree;
    }
    // todo: read the file listing paths to nested trees, 
    // call this constructor recursively on each, add to 
    // m_exclusive_trees
    m_vals = new GleamValues(vars, tuple);
}

double TreeFactory::Tree::operator()()const
{
    double test = (*m_filter_tree)(*m_vals);
    // either failed filter, or no nested trees
    if( m_exclusive_trees.size()==0 || test==0) return test;

    // now loop over the nested trees, return result of first non-zero evaluation
    for( std::vector<const TreeFactory::Tree*>::const_iterator it = m_exclusive_trees.begin();
        it!= m_exclusive_trees.end(); ++it)
    {
        test = (**it)();
        if( test>0) return test;
    }
    return 0; // no filter passed, or pure background.
}
std::string TreeFactory::Tree::title()const
{
    return m_filter_tree->title();
}
TreeFactory::Tree::~Tree()
{
    delete m_filter_tree;
    delete m_vals;
}
TreeFactory::~TreeFactory()
{
    ///@TODO: delete trees
}
