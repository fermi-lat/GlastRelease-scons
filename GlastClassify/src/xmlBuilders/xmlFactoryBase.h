/**@file xmlFactoryBase.h

@brief declaration of class xmFactoryBase
@author T. Usher

$Header$
*/

#ifndef xmlFactoryBase_h
#define xmlFactoryBase_h

#include "IxmlEngineFactory.h"
#include <xercesc/util/XercesDefs.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class  DOMElement;
XERCES_CPP_NAMESPACE_END
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

#include <vector>
#include <utility>
#include <map>
#include <iostream>

/** @class xmFactoryBase
@brief A factory for accessing decision trees

*/
class xmlFactoryBase 
{
public:
    // Useful typedefs...
    typedef std::vector<std::string> StringList;
    typedef std::vector<DOMElement*> DOMEvector;

    /** @brief constructotor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlFactoryBase(std::ostream& log=std::cout, int iVerbosity=0);

    DOMEvector  getXTPropertyVec(const DOMElement* element)                             const;
    DOMElement* getXTProperty(const DOMElement* element, const std::string& property)   const;
    DOMEvector  getXTSubPropVec(const DOMElement* element, const std::string& property) const;

    std::string getNextWord(std::string &sList, int &iEnd);
    
    double getNextDouble(std::string &sList, int &iEnd);

    std::string indent(int depth=0);

    virtual ~xmlFactoryBase();

    /// Finds first child given path
    const DOMElement* findXPath(const DOMElement* xmlParent, const std::vector<std::string>& nodeNames);

private:

    std::ostream& m_log;         //! output to this stream
    int           m_outputLevel; //! output level (verbosity)

};


#endif
