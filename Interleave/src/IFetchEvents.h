/** @file IFetchEvents.h

    @brief declaration of the IFetchEvents class

$Header$

*/

#ifndef IFetchEvents_h
#define IFetchEvents_h

#include <string>
class TChain;
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
    virtual int getFiles(double binVal, TChain* chain) = 0;

    /// Returns a TTree* constructed from the file stored in the m_dataStore
    virtual TTree* getTree(double binval){return 0;}

    virtual double minValFullRange()const{return -1e30;}///< return minimum value allowed
    virtual double maxValFullRange()const{return +1e30;}///< return maximum value allowed

    virtual double minVal()const{return -1e30; }///< return minimum value in current range
    virtual double maxVal()const{return +1e30; }///< return maximum value in current range

    /// test if value is valid in current range
    virtual bool isValid(double val)const{return val>=minVal() && val< maxVal();}

private:
    friend class XmlFetchEvents;

    /// path or file name associated with the dataStore which contains the information we wish to extract
    std::string m_dataStore; 
    /// name of the variable that we will use to search the dataStore
    std::string m_param;

};


#endif
