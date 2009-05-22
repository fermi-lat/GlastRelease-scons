/** @file ClassifyAlg.cxx
@brief Declaration and implementation of Gaudi algorithm ClassifyAlg

$Header$
*/

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include "classifier/DecisionTree.h"

#include "AtwoodTrees.h"
#include "GlastClassify/ITupleInterface.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T> class GleamItem : public GlastClassify::Item {
public:
    GleamItem<T>(const std::string& name, const std::string& type, T* data) :
      m_pdata(data), m_name(name), m_type(type)   {}
    
    virtual ~GleamItem<T>() {}

    operator double()const {
        return (double)*m_pdata;
    }

// LSR 14-Jul-08 code for ntuple types

    void setDataValue(void* data) 
    {
        if (m_type == "UInt_t")
        {
            *m_pdata = *(reinterpret_cast<int*>(data));
        }
        else if (m_type == "ULong64_t")
        {
            *m_pdata = *(reinterpret_cast<unsigned long long*>(data));
        }
        else if (m_type == "Float_t")
        {
            *m_pdata = *(reinterpret_cast<float*>(data));
        }
        else if (m_type == "Double_t")
        {
            *m_pdata = *(reinterpret_cast<double*>(data));
        }
        else if (m_type == "UChar_t")
        {
            memset(m_pdata, ' ', 80);
            strcpy(reinterpret_cast<char*>(m_pdata), reinterpret_cast<char*>(data));
        }
        else if (m_type == "Char_t")
        {
            memset(m_pdata, ' ', 80);
            strcpy(reinterpret_cast<char*>(m_pdata), reinterpret_cast<char*>(data));
        }
        else
        {
            throw std::invalid_argument("ClassifyAlg::Item: attempting to set an unrecognized data type");
        }
    }

private:
    T*          m_pdata;
    std::string m_name;
    std::string m_type;
    void*       m_treePtr;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class GleamTuple : public GlastClassify::ITupleInterface {
public:
    GleamTuple( INTupleWriterSvc* tuple, const std::string& treename)
        : m_tuple(tuple)
        , m_treename(treename)
    {
        m_tuple->getOutputTreePtr(m_treePtr, m_treename);

        if (!m_treePtr) 
            throw std::invalid_argument("GleamTuple constructor: can not get pointer to output tree!");

    }

// LSR 14-Jul-08 code for ntuple types

    const GlastClassify::Item* getItem(const std::string& name)const
    {
        const GlastClassify::Item* item = 0;
        void* dummy;

        std::string type = m_tuple->getItem(m_treename, name, dummy, m_treePtr);

        if (type == "Float_t")
        {
            float* data = (float*)dummy;
            item =new GleamItem<float>(name, type, data);
        }
        else if (type == "Double_t")
        {
            double* data = (double*)dummy;
            item =new GleamItem<double>(name, type, data);
        }
        else if (type == "UInt_t")
        {
            unsigned int* data = (unsigned int*)dummy;
            item =new GleamItem<unsigned int>(name, type, data);
        }
        else if (type == "ULong64_t")
        {
            unsigned long long* data = (unsigned long long*)dummy;
            item =new GleamItem<unsigned long long>(name, type, data);
        }
        else if (type == "Int_t")
        {
            int* data = (int*)dummy;
            item =new GleamItem<int>(name, type, data);
        }
        else if (type == "UChar_t")
        {
            char* data = (char*)dummy;
            item = new GleamItem<char>(name, type, data);
        }
        else if (type == "Char_t")
        {
            char* data = (char*)dummy;
            item = new GleamItem<char>(name, type, data);
        }
        else
        {
            throw std::invalid_argument("GleamTuple::getItem: attempting to set an unrecognized data type");
        }

        return item;
    }

// LSR 14-Jul-08 code for ntuple types

    /// create a new item (float only for now) in the tuple, which will take the given value
   void addItem(const std::string& name, float & value)
    {
        m_tuple->addItem(m_treename, name, &value);
    }
    
    void addItem(const std::string& name, double & value)
    {
        m_tuple->addItem(m_treename, name, &value);
    }
    void addItem(const std::string& name, unsigned long long & value)
    {
        m_tuple->addItem(m_treename, name, &value);
    }
    void addItem(const std::string& name, char & value)
    {
        m_tuple->addItem(m_treename, name, &value);
    }
private:
    INTupleWriterSvc*  m_tuple;
    void*              m_treePtr;
    const std::string& m_treename;
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class ClassifyAlg
@brief Extract info from tuple, etc. to add ft1 items to this of another tree
*/
class ClassifyAlg : public Algorithm {

public:
    ClassifyAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:
    StringProperty m_treename;
    StringProperty m_infoPath;
    StringProperty m_xmlFileName;

    /// this guy does the work!
    GlastClassify::AtwoodTrees * m_ctree;
    GleamTuple*                  m_tuple;
    bool                         m_treeInfo;

    int m_events;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static const AlgFactory<ClassifyAlg>  Factory;
const IAlgFactory& ClassifyAlgFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ClassifyAlg::ClassifyAlg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
,  m_ctree(0)
,  m_events(0)

{
    declareProperty("TreeName",      m_treename="MeritTuple");
    declareProperty("xmlFileName",   m_xmlFileName="");
    declareProperty("PrintTreeInfo", m_treeInfo=false);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode ClassifyAlg::initialize()
{
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // Use the Job options service to get the Algorithm's parameters
    setProperties();

    // get a pointer to RootTupleSvc 
    INTupleWriterSvc* rootTupleSvc(0);

    if( (sc = service("RootTupleSvc", rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
        return sc;
    }

    // create our interface to the tuple
    m_tuple= new GleamTuple(rootTupleSvc, m_treename);

    // create the classification object if requested
    try { 
        std::string path( m_xmlFileName.value()); 
        log << MSG::INFO;
        m_ctree = new  GlastClassify::AtwoodTrees(*m_tuple, log.stream(), path, m_treeInfo);
        log << endreq;
        
    }catch ( std::exception& e){
        log << MSG::ERROR << "Exception caught, class  "<< typeid( e ).name( ) << ", message:"
            << e.what() <<endreq;
        sc = StatusCode::FAILURE;
    }catch (...)  {
        log << MSG::ERROR << "Unexpected exception loading classification trees" << endreq;
        sc = StatusCode::FAILURE;
    }
    return sc;
}

StatusCode ClassifyAlg::execute() 
{
    MsgStream log(msgSvc(), name());
    if( m_ctree!=0){
        m_ctree->execute();
        ++m_events;
    }
    return StatusCode::SUCCESS;
}

StatusCode ClassifyAlg::finalize() 
{
    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "Processed " << m_events << " events." << endreq;
    delete m_ctree;
    delete m_tuple;
    setFinalized(); //  prevent being called again
    return StatusCode::SUCCESS;
}

