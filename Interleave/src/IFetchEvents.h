/** @file IFetchEvents.h

    @brief declaration of the IFetchEvents class

$Header$

*/

#ifndef IFetchEvents_h
#define IFetchEvents_h

#include <string>
#include "TChain.h"


/** @class IFetchEvents
    @brief manage the retrieval of events using some specified parameter(s)
    @author Heather Kelly heather625@gmail.com


*/
class IFetchEvents
{
public:

    /** @brief ctor
        @param param parameter used for retrieving data from xmlFile
        @param xmlFile name of the file containing the event store
    */
    IFetchEvents(const std::string& dataStore, const std::string& param) : m_dataStore(dataStore), m_param(param) {};

    ~IFetchEvents() { };

    /// Returns the value of the requested elemName, based on the bin associated with binVal
    virtual double getAttributeValue(const std::string& elemName, double binVal) = 0;

    virtual int getFiles(double binVal, TChain* chain) = 0;

private:
    friend class XmlFetchEvents;
    std::string m_dataStore; 
    std::string m_param; /// name of the variable that we use to key the rates

};


#endif
