/** @file RootTuple.cxx
    @brief implement class RootTuple

 $Header$
  Original author T. Burnett (w/ help from H. Kelley)
*/
#include "RootTuple.h"


// root includes
#include "TROOT.h"
#include "TFile.h"
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
#include <stdexcept>

namespace {
/*
    // convenient utility from Heather
TTree* getTree(TFile *f) {
  // Create an iterator on the list of keys
  TIter nextTopLevelKey(f->GetListOfKeys());
  TKey *keyTopLevel, *curDirKey;
  TTree* t=0; // return 

  // loop on keys, and search for the TTree named "t1"
  while  ( (keyTopLevel=(TKey*)nextTopLevelKey())!=0 ) {
    // I'm assuming we know the name of the TTree is "t1"
    TString name(keyTopLevel->GetName());
    TString className(keyTopLevel->GetClassName());

    if ((name.CompareTo("t1")==0) && (className.CompareTo("TTree")==0))  {
      // Found It
      t = (TTree*)f->Get(keyTopLevel->GetName());
      return t;
    }
    // If we find a directory - then we search it as well
    // Here I'm assuming that our directory structure only goes down one-level
    if (className.CompareTo("TDirectory")==0) {
      TDirectory *curDir = (TDirectory*)f->Get(name);
      TIter dirKeys(curDir->GetListOfKeys());
      while ( (curDirKey = (TKey*)dirKeys() ) ) {
        TString name(curDirKey->GetName());
        TString className(curDirKey->GetClassName());
        if ( (name.CompareTo("t1")==0) && (className.CompareTo("TTree")==0) ) {
          // Found it
          t = (TTree*)curDir->Get(curDirKey->GetName());
          return t;
        }
      }
    }
  }
  return t;
}
*/
class RootItem : public GlastClassify::Item {
public:
    RootItem(TLeaf * leaf) : m_leaf(leaf){}
    operator double() const { return m_leaf->GetValue();}
private:
    TLeaf* m_leaf;

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
    
    // Open the file, and get at the  TTree containing the data
    m_file =  new TFile(file.c_str(), "read");
    if( m_file==0 ) {
        std::cerr << "file \""<< file << "\" not found." << std::endl;
        throw std::runtime_error("File not found");
    }
    m_tree =  (TTree*)m_file->Get(treeName.c_str());
    if( m_tree ==0 ) {
        std::cerr << "tree \""<<treeName<< "\" not found." << std::endl;
        m_file->ls();
        throw std::invalid_argument("Tree not found");
    }
    m_numEvents = m_tree->GetEntries();
}

RootTuple::~RootTuple()
{
    if( m_output_file) m_output_file->Write();
    delete m_output_file;
}

const Item* RootTuple::getItem(const std::string& name)const
{
    TLeaf* leaf = m_tree->GetLeaf(name.c_str());
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
