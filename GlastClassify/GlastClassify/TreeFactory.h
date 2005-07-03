/**@file TreeFactory.h

@brief declaration of class TreeFactory

$Header$
*/

#ifndef GlastClassify_TreeFactory_h
#define GlastClassify_TreeFactory_h

class DecisionTree;


#include <string>
#include <vector>


namespace GlastClassify {

    /** @class TreeFactory
    @brief A factory for decision trees

    */
    class TreeFactory {
    public:

        /**
        nested class definition interface that must be implemented by client
        */
        class ILookupData {
        public:
            /// return pointer to data that will be filled by client later
            virtual const double * operator()(const std::string& name)=0;
            /// a bit of a kluge for ROOT. The above might have to be cast to a float*
            virtual bool isFloat()const {return false;}  
        };


        // forward declaration
        class GleamValues;

        TreeFactory(const std::string& path, ILookupData& lookup)
            :m_path(path)
            ,m_lookup(lookup){};

        /** @class Tree
        @brief nested class definition
        */
        class Tree {
        public:
            Tree( const std::string& path, ILookupData& lookup);
            double operator()()const;
            ~Tree();
            /// @return the title
            std::string title()const;

        private:
            const DecisionTree* m_dt;
            GleamValues* m_vals;
        };

        /// @return a reference to a new tree
        const TreeFactory::Tree& operator()(const std::string& name);


        /// @return value of Tree # i
        double evaluate(int i)const{return (*m_trees[i])();}

        ~TreeFactory();


    private:

        std::string m_path;
        ILookupData& m_lookup;

        std::vector<TreeFactory::Tree*> m_trees;

    };


} // namespace

#endif
