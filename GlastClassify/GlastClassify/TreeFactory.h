/**@file TreeFactory.h

@brief declaration of class TreeFactory
@author T. Burnett

$Header$
*/

#ifndef GlastClassify_TreeFactory_h
#define GlastClassify_TreeFactory_h

#if 0
#include "GlastClassify/ITreeFactory.h"
#endif
#include "GlastClassify/ITupleInterface.h"

class DecisionTree;

#include <vector>
#include <utility>


namespace GlastClassify {

    /** @class TreeFactory
    @brief A factory for accessing decision trees

    */
    class TreeFactory
#if 0 // THB, to isolate new structure for now
        : virtual public ITreeFactory 
#endif
    {
    public:

        // forward declaration
        class GleamValues;

        /** @brief ctor sets up for tree production 
        @param path file path to a folder containing tree data
        @param lookup Instance of a class supplied to by user, which is called back 
        to find address of each variable


        */
        TreeFactory(const std::string& path, ITupleInterface& tuple)
            :m_path(path)
            ,m_tuple(tuple){};

        /** @class Tree
        @brief nested class definition
        This class wraps a DecisionTree object
        */
        class Tree //: virtual public ITreeFactory::ITree  
        {
        public:
            Tree( const std::string& path, ITupleInterface& tuple);
            double operator()()const;
            ~Tree();
            /// @return the title
            std::string title()const;

        private:
            const DecisionTree* m_filter_tree; ///< optional filter, or filter and local tree
            std::vector<const TreeFactory::Tree*> m_exclusive_trees; ///< list of (assumed) exclusive trees to evaluate
            GleamValues* m_vals;
        };

        /** @param name a folder name completing the path to the folder containing the tree data
        
         @return a reference to a new tree. See also the evaluate() method.
         */
#if 0
        const ITreeFactory::ITree& operator()(const std::string& name);
#else
        const TreeFactory::Tree& operator()(const std::string& name);
#endif


        /// @return value of Tree # i for current set of values
        virtual double evaluate(int i) const; 

        /// index does the evaluate.
        virtual double operator[](int i) const{return evaluate(i);}

        virtual ~TreeFactory();


    private:

        std::string m_path;
        ITupleInterface& m_tuple;
#if 0
        std::vector<ITreeFactory::ITree*> m_trees;
#else
        std::vector<TreeFactory::Tree*> m_trees;
#endif

    };


} // namespace

#endif
