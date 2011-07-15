/** @file INTupleWriterSvc.h
    @brief declare abstract INtupleWriterSvc

    $Header$
*/
#ifndef _H_INTupleWriterSvc_
#define _H_INTupleWriterSvc_

#include "GaudiKernel/IInterface.h"
#include <string>

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_INTupleWriterSvc("INTupleWriterSvc",  9 ,0); 

/*! @class INTupleWriterSvc
 @brief Proper Gaudi abstract interface class for the ntupleWriterSvc 
*/
class INTupleWriterSvc : virtual public IInterface
{  

public:

    /// setup before event processing - required for Gaudi services
    virtual StatusCode initialize ()=0;
    
    /// cleanup after event processing - required for all Gaudi services
    virtual StatusCode finalize ()=0;

// LSR 14-Jul-08 code for ntuple types:

    /// add a pointer to a  double, or an array (depending on [n] following the name) of doubles
    virtual StatusCode addItem(const std::string & tupleName, 
             const std::string& itemName, const double* val,
             const std::string& fileName=std::string(""))=0;

    /// add a pointer to a  float, or an array of floats
    virtual StatusCode addItem(const std::string & tupleName, 
             const std::string& itemName, const float* val,
             const std::string& fileName=std::string(""))=0;
    
    /// add a pointer to an int (32 bits) or an array
    virtual StatusCode addItem(const std::string & tupleName, 
            const std::string& itemName, const int* val,
            const std::string& fileName=std::string(""))=0;
 
   
    /// add a pointer to a uint (32 bits) or an array
    virtual StatusCode addItem(const std::string & tupleName, 
            const std::string& itemName, const unsigned int* val,
            const std::string& fileName=std::string(""))=0;
 
    /// add a pointer to a uint (64 bits) or an array
    virtual StatusCode addItem(const std::string & tupleName, 
            const std::string& itemName, const unsigned long long* val,
            const std::string& fileName=std::string(""))=0;
 
        /** @brief Adds a pointer to a zero-terminated array of char 
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column.
    @param pval - pointer to the character array
     Note that the  data member must be the actual array containing the character data,
     which must be recopied for each event, say with strncpy.
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */

    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const char * pval,
            const std::string& fileName=std::string(""))=0;

#if 1 // deprecate! eliminate!

    virtual void storeRowFlag(bool flag)=0;
    virtual bool storeRowFlag()=0;
#endif
    /** store row flag by tuple Name option, retrive currrent
    @param tupleName Name of the tuple (TTree for RootTupleSvc implemetation)
    @param flag new value
    @return previous value
    If service does not implement, it is ignored (return false)
    */
    virtual bool storeRowFlag(const std::string& tupleName, bool flag)=0;

    /** @brief Set the pointer to the value of an existing item.
    @param tupleName  Name of the tuple
    @param itemName   Name of the item (or perhaps blank, see below)
    @param pointer    point pointer that will be set  

    Expect to throw exception if not found.
    @return a string describing the type (client then must cast)

    Note that the RootTupleSvc subclass implements a variation that it the name is empty,
    it will set the pointer to the TTree.
    
    */
    virtual std::string getItem(const std::string & tupleName, 
                                const std::string& itemName, void*& pointer,
                                const void* treePtr = 0) = 0;

    //! Provide access to input TTree pointer given tuple name
    virtual long long getInputTreePtr(void*& treePtr, 
                               const std::string& tupleName="MeritTuple") = 0;

    //! Provide access to output TTree pointer given tuple name
    virtual long long getOutputTreePtr(void*& treePtr, 
                              const std::string& tupleName="MeritTuple") = 0;

    //! Save the row in the output file
    virtual void saveRow(const std::string& tupleName)=0; 

    //! Returns merit version
    virtual int getMeritVersion() = 0;
    //! Set merit version
    virtual void setMeritVersion(int v) = 0;

    //! Allow clients to set basket size on a per branch basis, or for
    /// all branches by setting branchName="*" or branchName="xxx*"
   /// see TTree::SetBasketSize
    virtual void setBufferSize(const std::string& tupleName, int bufferSize, 
                         const std::string& branchName=std::string("*")) = 0;



    virtual bool getInputFileList(std::vector<std::string> &fileList) = 0;


    virtual bool setIndex( long long ) = 0 ;
    virtual long long index() = 0 ;
    virtual long long getNumberOfEvents() = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_INTupleWriterSvc; }

};


#endif // _H_INTupleWriterSvc
