/**@file TreeFactory.cxx

@brief implementation of class TreeFactory

$Header$
*/

#include "GlastClassify/TreeFactory.h"
#include "classifier/DecisionTree.h"
#include "classifier/TrainingInfo.h"
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
                throw std::invalid_argument("TreeFactory: did not find variable "+*it);
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

    std::string dtfile(path+"/dtree.txt");
    std::ifstream dtstream(dtfile.c_str());
    if(! dtstream.is_open()){
        throw std::runtime_error("failed to open file "+dtfile);
    }
    m_dt=new DecisionTree(dtstream);

    const TrainingInfo::StringList& vars = info.vars();
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
