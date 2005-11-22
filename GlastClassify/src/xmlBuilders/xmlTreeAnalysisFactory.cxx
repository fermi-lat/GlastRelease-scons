/**@file xmlTreeAnalysisFactory.cxx

@brief implementation of class xmlTreeAnalysisFactory

$Header$
*/

#include "xmlTreeAnalysisFactory.h"
#include "src/TreeAnalysis.h"

#include "src/XT/XTExprsnParser.h"

#include "xmlAppendEngineFactory.h"
#include "xmlCreateColumnsEngineFactory.h"
#include "xmlFilterColumnsEngineFactory.h"
#include "xmlFilterRowsEngineFactory.h"
#include "xmlMissingValuesEngineFactory.h"
#include "xmlModifyColumnsEngineFactory.h"
#include "xmlNewPredictEngineFactory.h"
#include "xmlReadTextFileEngineFactory.h"
#include "xmlShuffleEngineFactory.h"
#include "xmlSplitEngineFactory.h"
#include "xmlWriteTextFileEngineFactory.h"

#include <fstream>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <cmath>  // for M_PI, among others
    
/** @class Exception 
    @brief hold a string
*/
class Exception : public std::exception
{
public: 
    Exception(std::string error):m_what(error){}
    ~Exception() throw() {;}
    virtual const char *what( ) const  throw() { return m_what.c_str();} 
    std::string m_what;
};

#include "xmlBase/XmlParser.h"
#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"

XERCES_CPP_NAMESPACE_USE

typedef std::vector<DOMElement*> DOMEvector;

GlastClassify::xmlTreeAnalysisFactory::xmlTreeAnalysisFactory(const std::string& path, ITupleInterface& tuple)
                                                            : m_lookup(tuple), m_log(std::cout)
{
    //std::string sFileName = path+"/"+"DC2_Analysis_v2r1.imw";
    std::string sFileName = path;
    
    xmlBase::XmlParser xmlParser;

    m_domDocument = xmlParser.parse(sFileName.c_str());
    
    if(m_domDocument == 0)
    {
        // Error checking usually for a missing file.
        // Or sometimes due to a bad XML document.
        // Remember <P><BR></P> is only valid
        // in poorly formed html, not xml.
        // When we get a schema for UserLibrary, we'll
        // be able to use that as a validator also.
        std::invalid_argument("xmlTreeAnalysisFactory: could no open input file: " + sFileName);
    }

    return;
}

GlastClassify::TreeAnalysis* GlastClassify::xmlTreeAnalysisFactory::buildTreeAnalysis()
{
    GlastClassify::TreeAnalysis* tree = new GlastClassify::TreeAnalysis(m_lookup);

    // Clear everything just to be sure
    m_iNodeVec.clear();
    m_idToINodeMap.clear();
    m_typeToINodeVecMap.clear();

    // Find the activity nodes in the document
    int numNodes = findAllActivityNodes(tree);

    // Link them together
    int numLinks = linkActivityNodes(tree);

    return tree;
}

  
int GlastClassify::xmlTreeAnalysisFactory::findAllActivityNodes(GlastClassify::TreeAnalysis* tree)
{
    // Root...    
    DOMElement* domRoot = m_domDocument->getDocumentElement();

    // Build map between types and implementation of that type
    std::map<std::string, IxmlEngineFactory*> nodeFactoryMap;

    XTExprsnParser parser(tree->xtTupleMap());

    nodeFactoryMap["AppendEngineNode"]        = new xmlAppendEngineFactory(parser);
    nodeFactoryMap["CreateColumnsEngineNode"] = new xmlCreateColumnsEngineFactory(parser);
    nodeFactoryMap["FilterColumnsEngineNode"] = new xmlFilterColumnsEngineFactory(parser);
    nodeFactoryMap["FilterRowsEngineNode"]    = new xmlFilterRowsEngineFactory(parser);
    nodeFactoryMap["MissingValuesEngineNode"] = new xmlMissingValuesEngineFactory(parser);
    nodeFactoryMap["ModifyColumnsEngineNode"] = new xmlModifyColumnsEngineFactory(parser);
    nodeFactoryMap["PredictEngineNode"]       = new xmlNewPredictEngineFactory(parser);
    nodeFactoryMap["ReadTextFileEngineNode"]  = new xmlReadTextFileEngineFactory(parser);
    nodeFactoryMap["ShuffleEngineNode"]       = new xmlShuffleEngineFactory(parser);
    nodeFactoryMap["SplitEngineNode"]         = new xmlSplitEngineFactory(parser);
    nodeFactoryMap["WriteTextFileEngineNode"] = new xmlWriteTextFileEngineFactory(parser);

    // We now need to ensure that we're getting 
    // ActivityNode[@engineClass=='com.insightful.miner.PredictEngineNode']
    DOMEvector xmlActivityNodes;
    xmlBase::Dom::getDescendantsByTagName(domRoot, "ActivityNode", xmlActivityNodes);

    int cntem = 0;

    // Loop over the ActivityNodes in the input xml file
    for(DOMEvector::iterator actNodeIter = xmlActivityNodes.begin();
        actNodeIter != xmlActivityNodes.end(); actNodeIter++)
    {
        DOMElement* xmlActivityNode = *actNodeIter;

        std::string sType = xmlBase::Dom::getAttribute(xmlActivityNode, "engineClass");  
        std::string sId   = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

        // Get the label for this node
        DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
        std::string sDisplay    = xmlBase::Dom::getTagName(displayInfo);
        std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");

        // Parse the engine class a bit to get to the unique type identifier
        int iDelim = -1;
        std::string sNewType = getNextWord(sType, iDelim);

        while(iDelim > -1)
        {
            sNewType = getNextWord(sType, iDelim);
        }

        // Get this node
        IImActivityNode* iActivityNode = (*nodeFactoryMap[sNewType])(xmlActivityNode);

        cntem++;

        // Add this node to the list
        tree->addNewNode(iActivityNode);

        // Store in local vector and maps
        m_iNodeVec.push_back(iActivityNode);
        m_idToINodeMap[sId] = iActivityNode;
        m_typeToINodeVecMap[sNewType].push_back(iActivityNode);
    }

    //done
    return m_iNodeVec.size();
}


int GlastClassify::xmlTreeAnalysisFactory::linkActivityNodes(GlastClassify::TreeAnalysis* tree)
{
    // Root...    
    DOMElement* domRoot = m_domDocument->getDocumentElement();

    // We now need to ensure that we're getting 
    // ActivityNode[@engineClass=='com.insightful.miner.PredictEngineNode']
    std::vector<DOMElement *> xmlLinkVec;
    xmlBase::Dom::getDescendantsByTagName(domRoot, "Link", xmlLinkVec);

    int numLinks = xmlLinkVec.size();

    for(DOMEvector::iterator linkIter = xmlLinkVec.begin();
        linkIter != xmlLinkVec.end(); linkIter++)
    {
        DOMElement* xmlLink = *linkIter;

        std::string sFromNode = xmlBase::Dom::getAttribute(xmlLink, "fromNode");  
        std::string sFromPort = xmlBase::Dom::getAttribute(xmlLink, "fromPort");  
        int         fromPort  = xmlBase::Dom::getIntAttribute(xmlLink, "fromPort");  
        std::string sToNode   = xmlBase::Dom::getAttribute(xmlLink, "toNode");
        std::string sToPort   = xmlBase::Dom::getAttribute(xmlLink, "toPort");

        // Comment out for now... (fix in print statement?)
        //if (sFromPort == "1" && sToPort == "1") continue;

        IImActivityNode* fromINode = m_idToINodeMap[sFromNode];
        IImActivityNode* toINode   = m_idToINodeMap[sToNode];

        fromINode->setNodeLink(fromPort, toINode);

        // Find the "to" node in what remains of the vector of node pointers
        std::vector<IImActivityNode*>::iterator toNodeIter = std::find(m_iNodeVec.begin(), m_iNodeVec.end(), toINode);

        // Can have the case that a node has already been removed
        if (toNodeIter != m_iNodeVec.end())
        {
            m_iNodeVec.erase(toNodeIter);
        }
    }

    // Head vector should be last man standing?
    int headVecSize = m_iNodeVec.size();
    if (headVecSize > 0) 
    {
        IImActivityNode* headNode = m_iNodeVec.front();

        // Set the head node
        tree->setHeadNode(headNode);
    }
    else
    {
        throw Exception("ImSheetBuilder did not find a HEAD node!");
    }

    //done
    return numLinks;
}

// Helper function creates a std::string from a DOMString
// Usage:
// a = getNextWord(sList, iDelimPos=-1);
// b = getNextWord(sList, iDelimPos);

std::string GlastClassify::xmlTreeAnalysisFactory::getNextWord(std::string &sList, int &iEnd)
{
    // This allows space delimited data as well as comma, semicolon, etc.
    static char  cDELIM = '.'; 
    const char *cList;
    int iStart = iEnd + 1; // Start is previous end plus one.
    if(iStart > (int)sList.length())
    {
        // Serious error when you don't initialize iDelimPos=-1
        // before calling this function!
        throw Exception("Error with getNextWord!");
    }
    cList = sList.c_str();
    while((cList[iStart] == cDELIM) && cList[iStart] != 0)
    { 
        iStart++; 
    }
    iEnd = sList.find(cDELIM, iStart);
    if(iEnd > (int)sList.length())
    {
        throw Exception("Error with getNextWord");
    }
    return sList.substr(iStart, iEnd - iStart);
}
