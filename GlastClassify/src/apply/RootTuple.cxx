/** @file RootTuple.cxx
    @brief implement class RootTuple

 $Header$
  Original author T. Burnett (w/ help from H. Kelley)
*/
#include "RootTuple.h"


// root includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"
#include "TEventList.h"
#include "TTree.h"
#include "TSystem.h"
#include "TLeafF.h"
#include "TIterator.h"
#include "TKey.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TKey.h"
#include "TIterator.h"
#include "TString.h"

#include <cassert>
#include <stdexcept>
#include <iostream>

namespace {

class RootItem : public GlastClassify::Item {
public:
    RootItem(TLeaf * leaf) : m_leaf(leaf)
    {
        m_type  = m_leaf->GetTypeName();
        m_pdata = m_leaf->GetValuePointer();
    }
    operator double() const { return m_leaf->GetValue(); }

    void setDataValue(void* data) 
    {
        // Check that the data pointer is still valid (why do we need to do this?)
        void* tempPtr = m_leaf->GetValuePointer();

        // If it is not the same then we update.
// LSR 14-Jul-08 code for ntuple types
        if (m_pdata != tempPtr)
        {
            m_pdata = tempPtr;
        }

        if (m_type == "UInt_t")
        {
            *(reinterpret_cast<int*>(m_pdata)) = *(reinterpret_cast<int*>(data));
        }
        if (m_type == "ULong64_t")
        {
            *(reinterpret_cast<int*>(m_pdata)) = *(reinterpret_cast<unsigned long long*>(data));
        }
        else if (m_type == "Float_t")
        {
            *(reinterpret_cast<float*>(m_pdata)) = *(reinterpret_cast<float*>(data));
        }
        else if (m_type == "Double_t")
        {
            *(reinterpret_cast<double*>(m_pdata)) = *(reinterpret_cast<double*>(data));
        }
        else if (m_type == "UChar_t")
        {
            *(reinterpret_cast<std::string*>(m_pdata)) = 
                 *((reinterpret_cast<std::string*>(data))->data());
        }
        else
        {
            throw std::invalid_argument("ClassifyAlg::Item: attempting to set an unrecognized data type");
        }
    }
private:
    std::string m_type;
    void*       m_pdata;
    TLeaf*      m_leaf;
};

} // anonymous namespace
using namespace GlastClassify;

RootTuple::RootTuple( std::string file, std::string treeName)
:  m_event(0), m_output_tree(0), m_output_file(0) {

    // Initialize Root
    if ( 0 == gROOT )   {
        static TROOT meritRoot("root","ROOT I/O");
    } 
#ifdef WIN32
    int ret=gSystem->Load("libTree");
    if( ret==1) TTree dummy;

#endif
    // Use the following prescription for opening the input file 
    // in order to make compatible with xrootd
    //m_inChain = new TChain(treeName.c_str());

    //int rc = m_inChain->Add(file.c_str(), -1);
    //int en = m_inChain->LoadTree(0);
    //int rn = m_inChain->GetEntry(0);

    //int numEntries = m_inChain->GetEntries();

    // Open the file, and get at the  TTree containing the data
    //m_file =  new TFile(file.c_str(), "read");
    m_file =  TFile::Open(file.c_str(), "read");
    if( m_file==0 ) {
        std::cerr << "file \""<< file << "\" not found." << std::endl;
        throw std::runtime_error("File not found");
    }
    m_tree =  (TTree*)m_file->Get(treeName.c_str());
    //m_tree = m_inChain->GetTree();
    //m_tree = (TTree*)m_inChain;
    if( m_tree ==0 ) {
        std::cerr << "tree \""<<treeName<< "\" not found." << std::endl;
        //m_file->ls();
        throw std::invalid_argument("Tree not found");
    }
    m_numEvents = m_tree->GetEntries();

    // Prime the pump to set memory locations...
    m_tree->GetEvent(0);

    // Use this to define max tree size (apparantly a "bazillion" is not really a number...)
    Long64_t maxTreeSize=50000000000;

    TTree::SetMaxTreeSize(maxTreeSize);
}

RootTuple::~RootTuple()
{
    if( m_output_file) m_output_file->Write();
    delete m_output_file;
}

const Item* RootTuple::getItem(const std::string& name)const
{
    TLeaf* leaf = m_tree->GetLeaf(name.c_str());

    // Emulate nTupleWriterSvc by throwing an exception if we don't find 
    // a leaf for this variable
    if (leaf == 0)
    {
        throw std::invalid_argument("RootTuple::getItem cannot find leaf " + name);
    }
    return (leaf!=0) ? new RootItem(leaf) : 0;

}    

void RootTuple::addItem(const std::string& name, float& value)
{
    TLeaf* leaf = m_tree->GetLeaf(name.c_str());
    if( leaf!=0) {
        std::cout << "Adding item "<< name << ", which already exists" << std::endl;
        if( std::string(leaf->GetTypeName()) !="Float_t") {
            throw std::invalid_argument("RootTuple::addItem replacing wrong type");
        }
        leaf->SetAddress(&value);
    }else {
        m_tree->Branch(name.c_str(), (void*)&value, name.c_str());
    }
}
void RootTuple::addItem(const std::string& name, double& value)
{
    TLeaf* leaf = m_tree->GetLeaf(name.c_str());
    if( leaf!=0) {
        std::cout << "Adding item "<< name << ", which already exists" << std::endl;
        if( std::string(leaf->GetTypeName()) !="Double_t") {
            throw std::invalid_argument("RootTuple::addItem replacing wrong type");
        }
        leaf->SetAddress(&value);
    }else {
        m_tree->Branch(name.c_str(), (void*)&value, (name+"/D").c_str());
    }
}



bool RootTuple::nextEvent(){
    if(m_event<m_numEvents) {
        m_tree->GetEvent(m_event++);
        m_file->cd();
        return true;
    }
    return false;
}

    
bool RootTuple::getEvent(int idx)
{
    if (idx >= 0 && idx < m_numEvents)
    {
        m_tree->GetEvent(idx);
        return true;
    }
    return false;
}

void RootTuple::setOutputFile(const std::string & output_filename)
{
    m_output_file =new TFile(output_filename.c_str(), "recreate");

    m_output_tree = tree()->CloneTree(0); // do not copy events
}

void RootTuple::fill()
{
    assert(m_output_tree); 
    m_output_tree->Fill();
}
