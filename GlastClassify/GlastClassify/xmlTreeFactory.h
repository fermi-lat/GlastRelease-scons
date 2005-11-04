/**@file xmlTreeFactory.h

@brief declaration of class xmlTreeFactory
@author T. Burnett

$Header$
*/

#ifndef GlastClassify_xmlTreeFactory_h
#define GlastClassify_xmlTreeFactory_h

#include "GlastClassify/ITreeFactory.h"
#include <xercesc/util/XercesDefs.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class  DOMDocument;
XERCES_CPP_NAMESPACE_END
using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;

class DecisionTree;
//class DecisionTreeBuilder;
class ImSheetBuilder;

#include <vector>
#include <utility>
#include <map>

namespace GlastClassify {

    /** @class xmlTreeFactory
    @brief A factory for accessing decision trees

    */
    class xmlTreeFactory : virtual public ITreeFactory {
    public:

        // forward declaration
        class GleamValues;

        /** @brief ctor sets up for tree production 
        @param path file path to a folder containing tree data
        @param lookup Instance of a class supplied to by user, which is called back 
        to find address of each variable
        */
        xmlTreeFactory(const std::string& path, ITreeFactory::ILookupData& lookup);

        /** @param name a folder name completing the path to the folder containing the tree data
        
         @return a reference to a new tree. See also the evaluate() method.
         */
        const ITreeFactory::ITree& operator()(const std::string& name);


        /// @return value of Tree # i for current set of values
        virtual double evaluate(int i) const {return (*m_trees[i])();}

        /// index does the evaluate.
        virtual double operator[](int i) const{return evaluate(i);}

        virtual ~xmlTreeFactory();

    private:

        /** @class LocalTupleValues
        @brief This class calculates the local variables used in the xml Classification Trees
        */
        class LocalTupleValues
        {
        public:
            typedef std::map<std::string, double>        localValsMap;
            typedef std::map<std::string, const float*> nTupleVarsMap;

            LocalTupleValues(ITreeFactory::ILookupData& lookup);
            ~LocalTupleValues() {}

            bool isValue(const std::string& name) const 
                {return m_valsMap.find(name) != m_valsMap.end() ? true : false;}
            //const double* findVarPointer(std::string& name) const;

            double getValue(const std::string& name) const;
        private:
            localValsMap  m_valsMap;
            nTupleVarsMap m_tupleVals;
        };

        /** @class Tree
        @brief nested class definition
        This class wraps a DecisionTree object
        */
        class Tree : virtual public ITreeFactory::ITree {
        public:
            Tree(DecisionTree* tree, GleamValues* lookup) : m_dt(tree), m_vals(lookup) {}
            double operator()()const;
            ~Tree();
            /// @return the title
            std::string title()const;

        private:
            const DecisionTree* m_dt;
            GleamValues*        m_vals;
        };

        // The xml document containing the CT's
        DOMDocument*                       m_domDocument;

        // The "builder" for our Decision Trees
        //DecisionTreeBuilder*               m_builder;
        ImSheetBuilder*                    m_imSheet;

        // This looks up the values in the output ntuple
        ITreeFactory::ILookupData&         m_lookup;

        // Our collection of Classification Trees
        std::vector<xmlTreeFactory::Tree*> m_trees;

        // Class needed to calcluate local variables used in CT's
        LocalTupleValues                   m_localVals;

        // Mapping between Toby's CT naming convention and Bill's
        std::map<std::string,std::string>  m_TobyToBillMap;

        // Mapping between from IM's CT output variable names and Toby's
        std::map<std::string,std::string>  m_ImToTobyOutMap;

    };


} // namespace

#endif
