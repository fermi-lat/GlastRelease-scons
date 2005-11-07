/**@file IxmlEngineFactory.h

@brief declaration of class xmFactoryBase
@author T. Usher

$Header$
*/

#ifndef IxmlEngineFactory_h
#define IxmlEngineFactory_h

#include <xercesc/util/XercesDefs.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class  DOMElement;
XERCES_CPP_NAMESPACE_END
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

#include <vector>
#include <utility>
#include <map>
#include <iostream>

class IImActivityNode;

/** @class IxmlEngineFactory
@brief A factory 

*/
class IxmlEngineFactory 
{
public:
    /** @param name a folder name completing the path to the folder containing the tree data   
     @return a reference to a new tree. See also the evaluate() method.
     */
    virtual IImActivityNode* operator()(const DOMElement* element) = 0;
};


#endif
