/** @file IFetchEvents.h

    @brief declaration of the IFetchEvents class

$Header$

*/

#ifndef IFetchEvents_h
#define IFetchEvents_h

#include <string>
class TTree;


/** @class IFetchEvents
    @brief manage the retrieval of events using some specified parameter(s)
    @author Heather Kelly heather625@gmail.com


*/
class IFetchEvents
{
public:

    /** @brief ctor
        @param param parameter used for retrieving data from dataStore
        @param dataStore name of the file containing the event store
    */
    IFetchEvents(const std::string& dataStore, const std::string& param) : m_dataStore(dataStore), m_param(param) {};

    ~IFetchEvents() { };

    /// Returns the value of the requested elemName, based on the bin associated with binVal
    virtual double getAttributeValue(const std::string& elemName, double binVal) = 0;

    /// Returns a TChain* constructed from the fileList stored in the m_dataStore
    virtual int getFiles(double binVal, TTree* chain) = 0;

private:
    friend class XmlFetchEvents;

    /// path or file name associated with the dataStore which contains the information we wish to extract
    std::string m_dataStore; 
    /// name of the variable that we will use to search the dataStore
    std::string m_param;

};


#endif
