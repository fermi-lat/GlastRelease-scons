
#include "configData/base/ConfigTuple.h"

#include "TTree.h"
#include <iostream>

Int_t ConfigBranch::index(const ChannelKey& channel) const {
  Int_t idx = 0;
  Int_t factor = 1;
  if ( m_size.index3() > 1 ) {
    if ( channel.index3() >= 0 ) idx += channel.index3(); 
    factor *= m_size.index3();
  }
  if ( m_size.index2() > 1 ) { 
    if ( channel.index2() >= 0 ) idx += (channel.index2() * factor); 
    factor *= m_size.index2(); 
  }
  if ( m_size.index1() > 1 ) { 
    if ( channel.index1() >= 0 ) idx += (channel.index1() * factor); 
    factor *= m_size.index1(); 
  }
  if ( m_size.index0() > 1 ) { 
    if ( channel.index0() >= 0 ) idx += (channel.index0() * factor); 
  }
  return idx;
}

void ConfigBranch::setLeafSuffix(std::string& leafName) const {
  char tmp[8];
  if ( m_size.index0() > 1 ) { sprintf(tmp,"[%i]",m_size.index0()); leafName += tmp; }
  if ( m_size.index1() > 1 ) { sprintf(tmp,"[%i]",m_size.index1()); leafName += tmp; }
  if ( m_size.index2() > 1 ) { sprintf(tmp,"[%i]",m_size.index2()); leafName += tmp; }
  if ( m_size.index3() > 1 ) { sprintf(tmp,"[%i]",m_size.index3()); leafName += tmp; }
  leafName += "/"; leafName += m_type;
}

template <typename T>
void ConfigBranchImpl<T>::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  tree.SetBranchAddress(branchName.c_str(),(void*)(m_vals));
  if ( branchSize().singleton() ) return;
  std::string branchName_bcast = branchName;   branchName_bcast += "_bcast";
  tree.SetBranchAddress(branchName_bcast.c_str(),(void*)(&m_bcast));
}

template <typename T>
void ConfigBranchImpl<T>::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  std::string leafName = branchName;
  setLeafSuffix(leafName);
  tree.Branch(branchName.c_str(),(void*)(m_vals),leafName.c_str());
  if ( branchSize().singleton() ) return;
  std::string branchName_bcast = branchName;   branchName_bcast += "_bcast";
  std::string leafName_bcast = branchName_bcast;
  leafName_bcast += '/'; leafName_bcast += type();
  tree.Branch(branchName_bcast.c_str(),(void*)(&m_bcast),leafName_bcast.c_str());
  
}


void ConfigTuple::reset() {  
  for ( std::list<ConfigBranch*>::iterator itr = m_branchList.begin();
	itr != m_branchList.end(); itr++ ) {
    (*itr)->reset();
  }
}

void ConfigTuple::makeBranches(TTree& tree, const std::string& prefix) const {
  for ( std::list<ConfigBranch*>::const_iterator itr = m_branchList.begin();
	itr != m_branchList.end(); itr++ ) {
    (*itr)->makeBranch(tree,prefix);
  }
}

void ConfigTuple::attach(TTree& tree, const std::string& prefix) const {
  for ( std::list<ConfigBranch*>::const_iterator itr = m_branchList.begin();
	itr != m_branchList.end(); itr++ ) {
    (*itr)->attach(tree,prefix);
  }
}

void ConfigTuple::addBranch(ConfigBranch& aBranch) {
  std::map<std::string,ConfigBranch*>::iterator itrFind = m_branchMap.find(aBranch.name());
  if ( itrFind != m_branchMap.end() ) {
    std::cerr << "Branch " << aBranch.name() << " already exists" << std::endl;
    return;
  }
  m_branchMap[aBranch.name()] = &aBranch;
  m_branchList.push_back(&aBranch);
}

ConfigBranch* ConfigTuple::branch(const std::string& bname) {
  std::map<std::string,ConfigBranch*>::iterator itrFind = m_branchMap.find(bname);
  return itrFind != m_branchMap.end() ? itrFind->second: 0;
}

/// These are here so that CMT build the templates
/// Did I mention that I don't like CMT. 
ConfigBranchImpl<UShort_t> dummyShortBranch("",'s',ChannelKey());
ConfigBranchImpl<UInt_t> dummyIntBranch("",'i',ChannelKey());
ConfigBranchImpl<ULong64_t> dummyLongBranch("",'l',ChannelKey());

