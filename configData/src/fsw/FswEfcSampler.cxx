//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      A GEM ROI register class 
//
// Environment:
//      Software developed for GLAST.
//
// Author List:
//      Martin Kocian
//
// Copyright Information:
//      Copyright (C) 2005      Stanford Linear Accelerator Center
//
//---------------------------------------------------------------------------

#include "configData/fsw/FswEfcSampler.h"
#include <iomanip>
#include <algorithm> 

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include "facilities/Util.h"


FswEfcSampler* FswEfcSampler::makeFromXmlFile(const char* fileName) {
  static xmlBase::XmlParser parser(true);
  xmlBase::DOMDocument* doc(0);
  try {
    doc =parser.parse( fileName );
  }
  catch (xmlBase::ParseException ex) {
    std::cerr << "caught exception with message " << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    return 0;
  }
  xmlBase::DOMElement* elt = doc->getDocumentElement();
  if ( elt == 0 ) return 0;

  std::vector<xmlBase::DOMElement*> enabledList;
  std::vector<xmlBase::DOMElement*> inputPSList;
  std::vector<xmlBase::DOMElement*> outputPSList;
  std::vector<xmlBase::DOMElement*> psList;
  xmlBase::Dom::getDescendantsByTagName(elt,"Enabled",enabledList);
  xmlBase::Dom::getDescendantsByTagName(elt,"Input_Prescale",inputPSList);
  xmlBase::Dom::getDescendantsByTagName(elt,"Output_Prescale",outputPSList);
  xmlBase::Dom::getDescendantsByTagName(elt,"Prescale",psList);
  if ( enabledList.size() != 1 ||
       inputPSList.size() != 1 ||
       outputPSList.size() != 1 ||
       psList.size() != 32 ) return 0;  
  std::string enabledSt = xmlBase::Dom::getTextContent( enabledList[0] );  
  std::string inputSt = xmlBase::Dom::getTextContent( inputPSList[0] );
  std::string outputSt = xmlBase::Dom::getTextContent( outputPSList[0] );

  FswEfcSampler* retVal(0);
  try {
    unsigned enabled = facilities::Util::stringToUnsigned( enabledSt );
    unsigned input = facilities::Util::stringToUnsigned( inputSt );
    unsigned output = facilities::Util::stringToUnsigned( outputSt );
    unsigned ps[32];    
    for ( int i(0); i < 32; i++ ) {
      std::string psSt = xmlBase::Dom::getTextContent(psList[i]);
      int lineNum = xmlBase::Dom::getIntAttribute(psList[i],"line");
      ps[lineNum] = facilities::Util::stringToUnsigned(psSt);
    }
    retVal = new FswEfcSampler;
    retVal->set(ps,input,output,enabled);
  } catch (...) {
    delete retVal;
    return 0;
  }
  return retVal;
}

unsigned FswEfcSampler::mapKey(unsigned lpaMode, unsigned handlerId) {
  return (lpaMode << 8) | handlerId;
}

FswEfcSampler::FswEfcSampler()
  :ConfigBranch("Prescales",'i',ChannelKey(32)){  
  clear();
}

void FswEfcSampler::clear(){
  for (int i=0;i<32;i++){
    m_prescalers[i]=0;
  }
  m_inputPrescaler = 0;
  m_outputPrescaler = 0;
  m_enabled = 0;
}

void FswEfcSampler::set(unsigned prescalers[32], unsigned input, unsigned output, unsigned enabled){
  for ( int i(0); i < 32; i++ ) {
    m_prescalers[i] = prescalers[i];    
  }
  m_inputPrescaler = input;
  m_outputPrescaler = output;
  m_enabled = enabled;
}

unsigned FswEfcSampler::prescaleFactor(enums::Lsf::RsdState rsdState, enums::Lsf::LeakedPrescaler leakedPrescaler) const {
  switch ( rsdState ) {
  case enums::Lsf::VETOED:
    return 0;
  case enums::Lsf::PASSED:
    return 1;
  case enums::Lsf::IGNORED:
    return m_inputPrescaler;
  case enums::Lsf::INVALID:
    return LSF_INVALID_UINT;
  case enums::Lsf::SUPPRESSED:
  case enums::Lsf::LEAKED:
    if ( leakedPrescaler >= 0 ) { return m_prescalers[leakedPrescaler]; }
    else if ( leakedPrescaler == enums::Lsf::OUTPUT ) { return m_outputPrescaler; }
    else if ( leakedPrescaler == enums::Lsf::INPUT ) { return m_inputPrescaler; }
    break;
  }
  return LSF_INVALID_UINT;  
}


// Attach this value to a TTree
void FswEfcSampler::makeBranch(TTree& /* tree */, const std::string& /* prefix */) const {
  return;
}

void FswEfcSampler::attach(TTree& tree, const std::string& prefix) const {
  return;
}

std::ostream& operator <<(std::ostream& os,  const FswEfcSampler& tc){
  os<<"Input Prescale: " << tc.inputPrescaler() <<std::endl;
  os<<"Output Prescale: " << tc.outputPrescaler() <<std::endl;
  os<<"Prescales: "<<std::endl;
  for (int i=0;i<32;i++){
    os<< "  line " << std::setw(2)<< i <<": " << tc.prescaler(i);
    os<< ( (tc.enabled() & (1 << i)) ?  " On" : " Off") << std::endl;
  }
  return os;
}

