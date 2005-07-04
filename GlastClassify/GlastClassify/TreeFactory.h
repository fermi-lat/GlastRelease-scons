/**@file TreeFactory.h

@brief declaration of class TreeFactory
@author T. Burnett

$Header$
*/

#ifndef GlastClassify_TreeFactory_h
#define GlastClassify_TreeFactory_h

class DecisionTree;


#include <string>
#include <vector>


namespace GlastClassify {

    /** @class TreeFactory
    @brief A factory for accessing decision trees

    */
    class TreeFactory {
    public:

        /** @class ILookupData
        @brief nested class definition interface that must be implemented by client
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

        /** @brief ctor sets up for tree production 
        @param path file path to a folder containing tree data
        @param lookup Instance of a class supplied to by user, which is called back 
        to find address of each variable


        */
        TreeFactory(const std::string& path, ILookupData& lookup)
            :m_path(path)
            ,m_lookup(lookup){};

        /** @class Tree
        @brief nested class definition
        This class wraps a DecisionTree object
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

        /** @param name a folder name completing the path to the folder containing the tree data
        
         @return a reference to a new tree. See also the evaluate() method.
         */
        const TreeFactory::Tree& operator()(const std::string& name);


        /// @return value of Tree # i for current set of values
        double evaluate(int i)const{return (*m_trees[i])();}

        ~TreeFactory();


    private:

        std::string m_path;
        ILookupData& m_lookup;

        std::vector<TreeFactory::Tree*> m_trees;

    };


} // namespace

#endif
