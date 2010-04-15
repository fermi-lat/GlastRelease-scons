
#include "facilities/Util.h"
#include <sstream>

#include "configData/ConfigReader.h"

#include "configData/base/ConfigTuple.h"
#include "configData/gem/TrgConfigParser.h"
#include "configData/db/MootDBImpl.h"
#include "configData/db/LatcDBImpl.h"

#include "commonRootData/AcdMap.h"
#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"
#include <xercesc/dom/DOMNode.hpp>

using namespace xmlBase;
#ifdef WIN32 // this is ugly, apparently gcc has such a define.
# define or ||
#endif

namespace{
    char ** nullPtr = 0;
}

// local function to return the bcast and non-bcast LATC fragments, given a MOOT key
bool _latcFromMootKey(const unsigned configKey, 
                      std::vector<std::string>& bcastFiles,
                      std::vector<std::string>& latcFiles){
  std::vector<std::string> filelist;
  MootDBImpl moot;
  filelist=moot.getFilenameList(configKey);
  for ( std::vector<std::string>::iterator it=filelist.begin();it!=filelist.end();it++){
    std::string latcFileFull=*it;
    if ( latcFileFull.find("DFT") !=  latcFileFull.npos ) {
      bcastFiles.push_back(latcFileFull);
    } else {
      latcFiles.push_back(latcFileFull);
    }
  }
  return kTRUE;
}

// local function to return the bcast and non-bcast LATC fragments, given a LATC master key
bool _latcFromMaster(const unsigned masterKey, 
                      std::vector<std::string>& bcastFiles,
                      std::vector<std::string>& latcFiles){
  std::vector<std::string> filelist;
  LatcDBImpl latc;
  filelist=latc.getFilenameList(masterKey);
  for ( std::vector<std::string>::iterator it=filelist.begin();it!=filelist.end();it++){
    std::string latcFileFull=*it;
    if ( latcFileFull.find("DFT") !=  latcFileFull.npos ) {
      bcastFiles.push_back(latcFileFull);
    } else {
      latcFiles.push_back(latcFileFull);
    }
  }
  return kTRUE;
}

// local function to return bcast and non-bcast LATC frags given an intSeApp style app XML
bool _latcFromAppXML(const std::string& fileName,
                     const std::string& latcPath,
                     std::vector<std::string>& bcastFiles,
                     std::vector<std::string>& latcFiles){
  XmlParser parser( true );
  DOMDocument* doc(0);
  try {
    doc = parser.parse( fileName.c_str() );
  } 
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    return kFALSE;
  }
    
  DOMElement* topElt = doc->getDocumentElement();
  std::vector<DOMElement*> eltList;
    
    
  std::vector<DOMElement*> everythingList;
    
  Dom::getDescendantsByTagName( topElt, "everything", everythingList );
  for ( std::vector<DOMElement*>::iterator itrEv = everythingList.begin();
        itrEv != everythingList.end(); ++itrEv ) {
    std::vector<std::string> fileNames;
    std::string everything = Dom::getTextContent(*itrEv);    
    facilities::Util::stringTokenize(everything,",",fileNames);
    for ( std::vector<std::string>::const_iterator itrFileName = fileNames.begin(); itrFileName != fileNames.end(); itrFileName++ ) {
      std::string latcFileFull = latcPath; latcFileFull += '/'; latcFileFull += (*itrFileName);
      if ( latcFileFull.find("DFT") !=  latcFileFull.npos ) {
        bcastFiles.push_back(latcFileFull);
      } else {
        latcFiles.push_back(latcFileFull);
      }
    }
  }
    
  if ( everythingList.size() == 0 ) {
    Dom::getDescendantsByTagName( topElt, "paramFile", eltList );
  }
    
  std::vector<DOMElement*> latcNodeList;
  Dom::getDescendantsByTagName( topElt, "latcFiles", latcNodeList);
    
  for ( std::vector<DOMElement*>::iterator itrLatc = latcNodeList.begin();
        itrLatc != latcNodeList.end(); ++itrLatc ) {
    Dom::getDescendantsByTagName(*itrLatc,"latcFile", eltList , false);
    Dom::getDescendantsByTagName(*itrLatc,"TFE",eltList,false);
    Dom::getDescendantsByTagName(*itrLatc,"TDC",eltList,false);
    Dom::getDescendantsByTagName(*itrLatc,"TEM",eltList,false);
    Dom::getDescendantsByTagName(*itrLatc,"CFE",eltList,false);
    Dom::getDescendantsByTagName(*itrLatc,"AFE",eltList,false);
    Dom::getDescendantsByTagName(*itrLatc,"ARC",eltList,false);
    Dom::getDescendantsByTagName(*itrLatc,"SPT",eltList,false);
    Dom::getDescendantsByTagName(*itrLatc,"bcast",eltList,false);    
  }
    
  std::vector<DOMElement*>::iterator itr = eltList.begin();  

  for ( itr = eltList.begin(); itr != eltList.end(); itr++ ) {      
    std::string latcFileBase = Dom::getTextContent(*itr);
    std::string latcFileFull = latcPath; latcFileFull += '/'; latcFileFull += latcFileBase;
    if ( latcFileBase.find("DFT") !=  latcFileBase.npos ) {
      bcastFiles.push_back(latcFileFull);
    } else {
      latcFiles.push_back(latcFileFull);
    }
  }
  return kTRUE;
}

Bool_t ConfigReader::readTopLvl(const std::string& fileName, 
                                const std::string& latcPath, 
                                bool bcast, bool useMootKey) {
  std::vector<std::string> bcastFiles;
  std::vector<std::string> latcFiles;
  char *endptr;
  unsigned key=strtoul(fileName.c_str(),&endptr,0);
  if (endptr[0]=='\0'){ // This is a database key
    if (useMootKey) {
      if (! _latcFromMootKey(key, bcastFiles, latcFiles))
        return kFALSE;
    } else {
      if (! _latcFromMaster(key, bcastFiles, latcFiles))
        return kFALSE;
    }
  }else{
    if (! _latcFromAppXML(fileName, latcPath, bcastFiles, latcFiles))
      return kFALSE;
  }
    
  std::vector<std::string>::const_iterator itrF = bcastFiles.begin();
  for (  itrF = bcastFiles.begin(); itrF != bcastFiles.end(); itrF++ ) {
    if ( ! read(*itrF) ) {
      std::cout << "Failed to read " << *itrF << std::endl;
    }else{
      //      std::cout << "Read bcast file " << *itrF << std::endl;
    }
  }
  if ( !bcast ) {
    for (  itrF = latcFiles.begin(); itrF != latcFiles.end(); itrF++ ) {
      if ( ! read(*itrF) ) {
	std::cout << "Failed to read " << *itrF << std::endl;
      }else{
	      //std::cout << "Read file" << *itrF << std::endl;
      }
    }
  }
  loadAcdIds();

  return kTRUE;

}

Bool_t ConfigReader::read(const std::string& fileName) {
  XmlParser parser( true );
  DOMDocument* doc(0);
  try {
    doc = parser.parse( fileName.c_str() );
  } 
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    return kFALSE;
  }
  DOMElement* topElt = doc->getDocumentElement();
  
  std::vector<DOMElement*> eltList;

  Dom::getChildrenByTagName( topElt, "AEM", eltList );
  std::vector<DOMElement*>::iterator itr = eltList.begin();  

  for ( itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( read_AEM(**itr) == kFALSE ) return kFALSE;    
  }
  Dom::getChildrenByTagName( topElt, "TEM", eltList );
  for ( itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( read_TEM(**itr) == kFALSE ) return kFALSE;    
  }
  Dom::getChildrenByTagName( topElt, "GEM", eltList );
  for ( itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( read_GEM(**itr) == kFALSE ) return kFALSE;    
  }
  return kTRUE;
}

Bool_t ConfigReader::read_AEM(DOMElement& elem) {
  const std::string ARC("ARC");  
  ChannelKey key(0);
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,ARC) ) {
      if ( read_ARC(**itr) == kFALSE ) return kFALSE;    
    } else {
      if ( ! readUInt(**itr,key) ) {
	continue;
      }
    }
  }
  return kTRUE;
}

Bool_t ConfigReader::read_ARC(DOMElement& elem) {  
  const std::string AFE("AFE");
  const std::string pha("pha_threshold");
  if ( ! ConfigReader::getId(elem,_iARC) ) return kFALSE;
  ChannelKey key(_iARC);
  std::vector<DOMElement*> eltList;
  //std::cout << Dom::getTagName(&elem) << std::endl;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,AFE) ) {
      if ( read_AFE(**itr) == kFALSE ) return kFALSE;
      continue;
    }  else if ( Dom::getTagName(*itr).find(pha) == ! std::string::npos ) {
      Int_t idx = getIndex(**itr,pha);
      ChannelKey keyPha(_iARC,idx);
      if ( ! readUShort_pha(**itr,keyPha,pha.c_str()) ) {
	continue;
      }
    } else {
      if ( ! readUShort(**itr,key) ) {
	continue;
      }
    } 
  } 
  return kTRUE;
}

Bool_t ConfigReader::read_AFE(DOMElement& elem) {
  if ( ! ConfigReader::getId(elem,_iAFE) ) return kFALSE;
  ChannelKey key(_iARC,_iAFE);  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  //std::cout << Dom::getTagName(&elem) << std::endl;
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( ! readUShort(**itr,key) ) {
      continue;
    }
  }
  return kTRUE;
}


Bool_t ConfigReader::read_TEM(DOMElement& elem){
  if ( ! ConfigReader::getId(elem,_iTEM) ) return kFALSE;
  ChannelKey key(_iTEM);  
  const std::string SPT("SPT");
  const std::string TIC("TIC");
  const std::string TCC("TCC");
  const std::string CCC("CCC");  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,SPT) ) {
      if ( read_SPT(**itr) == kFALSE ) return kFALSE;    
    } else if  ( Dom::checkTagName(*itr,TIC) ) {
      if ( read_TIC(**itr) == kFALSE ) return kFALSE; 
    } else if  ( Dom::checkTagName(*itr,TCC) ) {
      if ( read_TCC(**itr) == kFALSE ) return kFALSE; 
    } else if  ( Dom::checkTagName(*itr,CCC) ) {
      if ( read_CCC(**itr) == kFALSE ) return kFALSE; 
    } else {
      if ( ! readUInt(**itr,key) ) {
	continue;
      }
    }
  }
  return kTRUE;
}

Bool_t ConfigReader::read_TIC(DOMElement& elem){
  ChannelKey key(_iTEM);
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( ! readUInt(**itr,key) ) {
      continue;
    }
  }
  return kTRUE;
}


Bool_t ConfigReader::read_CCC(DOMElement& elem){
  if ( ! ConfigReader::getId(elem,_iCCC) ) return kFALSE;
  ChannelKey key(_iTEM,_iCCC);  
  const std::string CRC("CRC");  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,CRC) ) {
      if ( read_CRC(**itr) == kFALSE ) return kFALSE;    
    } else {
      if ( ! readUInt(**itr,key) ) {
	continue;
      }
    }
  }
  return kTRUE;
};

Bool_t ConfigReader::read_CRC(DOMElement& elem){
  if ( ! ConfigReader::getId(elem,_iCRC) ) return kFALSE;
  ChannelKey key(_iTEM,_iCCC,_iCRC);  
  const std::string CFE("CFE");  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,CFE) ) {
      if ( read_CFE(**itr) == kFALSE ) return kFALSE;    
    } else {
      if ( ! readUInt(**itr,key) ) {
	continue;
      }
    }
  }
  return kTRUE;
};

Bool_t ConfigReader::read_CFE(DOMElement& elem){
  if ( ! ConfigReader::getId(elem,_iCFE) ) return kFALSE;
  ChannelKey key(_iTEM,_iCCC,_iCRC,_iCFE);  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  //std::cout << Dom::getTagName(&elem) << std::endl;
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( ! readUInt(**itr,key) ) {
      continue;
    }
  }
  return kTRUE;
};

Bool_t ConfigReader::read_SPT(DOMElement& elem){
  if ( ! ConfigReader::getSptId(elem,_iSPT) ) return kFALSE;
  ChannelKey key(_iTEM,_iSPT);  
  const std::string TFE("TFE");  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,TFE) ) {
      if ( read_TFE(**itr) == kFALSE ) return kFALSE;    
    } else {
      if ( ! readUInt(**itr,key) ) {
	continue;
      }
    }
  }
  return kTRUE;
};

Bool_t ConfigReader::read_TFE(DOMElement& elem){
  std::string spt = Dom::getAttribute(&elem,"ID");
  if(spt.find("x") == spt.npos && spt.find("y") == spt.npos){
    if ( ! ConfigReader::getId(elem,_iTFE) ) {
      std::cerr << "couldn't find ID for tfe " << _iTFE << std::endl;
      return kFALSE;
    }
  }else{
    if ( ! ConfigReader::getTfeId(elem,_iTFE) ) {
      std::cerr << "couldn't find text ID for tfe " << _iTFE << std::endl;
      return kFALSE;
    }
  }
  ChannelKey key(_iTEM,_iSPT,_iTFE);  
  const std::string TDC("TDC");  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  //std::cout << Dom::getTagName(&elem) << std::endl;
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,TDC) ) {
      if ( read_TDC(**itr) == kFALSE ) {
	std::cerr << "failed to read tdc " << _iTFE << std::endl;
	return kFALSE;    
      }
    } else {
      if ( ! readULong(**itr,key) ) {
	std::cerr << "failed to long for " <<  _iTFE << std::endl;
	continue;
      }
    }
  }
  return kTRUE;
};

Bool_t ConfigReader::read_TDC(DOMElement& elem){
  ChannelKey key(_iTEM,_iSPT,_iTFE);
  DOMElement* tfe_dac = Dom::findFirstChildByName( &elem, "tfe_dac");
  std::vector<DOMElement*> eltList;  
  //std::cout << Dom::getTagName(&elem) << std::endl;  
  if ( tfe_dac == 0 ) return kTRUE;
  Dom::getChildrenByTagName( tfe_dac, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( ! readUInt(**itr,key) ) {
      continue;
    }
  }
  return kTRUE;
};

Bool_t ConfigReader::read_TCC(DOMElement& elem){
  if ( ! ConfigReader::getId(elem,_iTCC) ) return kFALSE;
  ChannelKey key(_iTEM,_iTCC);  
  const std::string TRC("TRC");  
  std::vector<DOMElement*> eltList;
  //std::cout << Dom::getTagName(&elem) << std::endl;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( Dom::checkTagName(*itr,TRC) ) {
      if ( read_TRC(**itr) == kFALSE ) return kFALSE;    
    } else {
      if ( ! readUInt(**itr,key) ) {
	continue;
      }
    }
  }
  return kTRUE;
    
};

Bool_t ConfigReader::read_TRC(DOMElement& elem){
  if ( ! ConfigReader::getId(elem,_iTRC) ) return kFALSE;
  ChannelKey key(_iTEM,_iTCC,_iTRC);  
  std::vector<DOMElement*> eltList;
  Dom::getChildrenByTagName( &elem, "*", eltList );
  //std::cout << Dom::getTagName(&elem) << std::endl;
  for ( std::vector<DOMElement*>::iterator itr = eltList.begin(); itr != eltList.end(); itr++ ) {
    if ( ! readULong(**itr,key) ) {
      continue;
    }
  }
  return kTRUE;
};

Bool_t ConfigReader::read_GEM(DOMElement& elem){
  TrgConfigParser parser(&elem);
  parser.parse(m_gem);
  return kTRUE;
};

Bool_t ConfigReader::getId(DOMElement& elem, Int_t& id) {
  id = ChannelKey::UNDEF;
  std::string ID("ID");
  std::string BCAST("BCAST");
  std::string bcast("bcast");
  if ( Dom::getAttribute(&elem,ID) == BCAST || Dom::getAttribute(&elem,ID) == bcast ) {
    id = ChannelKey::BCAST;
    return kTRUE;    
  }
  id = Dom::getIntAttribute(&elem,ID);
  return kTRUE;
}

Bool_t ConfigReader::getUShort(DOMElement& elem, UShort_t& val) {
  std::vector<DOMElement*> rList;
  Dom::getChildrenByTagName( &elem, "register", rList );
  std::string str;
  if(rList.empty()){  
    str = Dom::getTextContent( &elem );
  }else{
    str=Dom::getTextContent(rList[0]);
  }
  // should add check
  val = (UShort_t)(strtoul(str.c_str(),nullPtr,0));
  return kTRUE;
}

Bool_t ConfigReader::getUInt(DOMElement& elem, UInt_t& val) {
  std::vector<DOMElement*> rList;
  Dom::getChildrenByTagName( &elem, "register", rList );
  std::string str;
  if(rList.empty()){  
    str= Dom::getTextContent(&elem);
  }else{
    str=Dom::getTextContent(rList[0]);
  }
  // should add check
  val = (UInt_t)strtoul(str.c_str(),nullPtr,0);
  return kTRUE;
}

Bool_t ConfigReader::getULong(DOMElement& elem, ULong64_t& val) {
  std::vector<DOMElement*> rList;
  Dom::getChildrenByTagName( &elem, "register", rList );
  std::string str;
  if(rList.empty()){  
    str = Dom::getTextContent(&elem);
  }else{
    str=Dom::getTextContent(rList[0]);
  }
  // should add check
//  val = (ULong64_t)strtoull(str.c_str(),nullPtr,0);
  val = (ULong64_t)facilities::Util::stringToUll(str);
  return kTRUE;
}


Bool_t ConfigReader::readUShort(DOMElement& elem, const ChannelKey& key, const char* bName) {  
  const std::string tagName = bName == 0 ? Dom::getTagName(&elem) : bName;
  ConfigBranchImpl<UShort_t>* reg = static_cast<ConfigBranchImpl<UShort_t>*>(m_config->branch(tagName));
  UShort_t val(0);
  if ( reg == 0 || ! getUShort(elem,val) ) {
    return kFALSE;
  }
  if (key.hasBCAST()) {
    reg->setAll(val);
  } else {
    reg->setVal(key,val);
  }
  return kTRUE;
}

Bool_t ConfigReader::readUShort_pha(DOMElement& elem, const ChannelKey& key, const char* bName) {  
  const std::string tagName = bName == 0 ? Dom::getTagName(&elem) : bName;
  ConfigBranchImpl<UShort_t>* reg = static_cast<ConfigBranchImpl<UShort_t>*>(m_config->branch(tagName));
  UShort_t val(0);
  if ( reg == 0 || ! getUShort(elem,val) ) {
    return kFALSE;
  }
  if (key.hasBCAST()) {
    for ( UInt_t garc(0); garc < 12; garc++ ) {
      ChannelKey bcastkey(garc,key.index1());
      reg->setVal(bcastkey,val);
    }
  } else {
    reg->setVal(key,val);
  }
  return kTRUE;
}

Bool_t ConfigReader::readUInt(DOMElement& elem, const ChannelKey& key, const char* bName ) {
  const std::string tagName =  bName == 0 ? Dom::getTagName(&elem) : bName;
  ConfigBranchImpl<UInt_t>* reg = static_cast<ConfigBranchImpl<UInt_t>*>(m_config->branch(tagName));
  UInt_t val(0);
  if ( reg == 0 || ! getUInt(elem,val) ) {
    return kFALSE;
  }
  if (key.hasBCAST()) {
    reg->setAll(val);
  } else {
    reg->setVal(key,val);
  }
  return kTRUE;
}

Bool_t ConfigReader::readULong(DOMElement& elem, const ChannelKey& key, const char* bName ) {
  const std::string tagName =  bName == 0 ? Dom::getTagName(&elem) : bName;
  ConfigBranchImpl<ULong64_t>* reg = static_cast<ConfigBranchImpl<ULong64_t>*>(m_config->branch(tagName));
  ULong64_t val(0);
  if ( reg == 0 || ! getULong(elem,val) ) {
    return kFALSE;
  }
  if (key.hasBCAST()) {
    reg->setAll(val);
  } else {
    reg->setVal(key,val);
  }
  return kTRUE;
}

Int_t ConfigReader::getIndex(DOMElement& elem, const std::string& rootName) {
  std::string fullName = Dom::getTagName(&elem);
  std::string endName(fullName,rootName.length()+1);
  return strtol(endName.c_str(),nullPtr,10);
}

Bool_t ConfigReader::getSptId(DOMElement& elem, Int_t& id) {
  id = ChannelKey::UNDEF;
  std::string ID("ID");
  std::string spt = Dom::getAttribute(&elem,ID);
  std::string BCAST("BCAST");
  std::string bcast("bcast");
  if ( spt == BCAST || spt == bcast) { 
    id = ChannelKey::BCAST;
    return kTRUE;    
  }
  static std::map<std::string,int> sptMap;
  if ( sptMap.size() == 0 ) {
    sptMap["+x0"] = 0; sptMap["+y0"] = 9; sptMap["-y0"] = 18; sptMap["-x0"] = 27;
    sptMap["+x1"] = 1; sptMap["+y1"] = 10; sptMap["-y1"] = 19; sptMap["-x1"] = 28;
    sptMap["+x2"] = 2; sptMap["+y2"] = 11; sptMap["-y2"] = 20; sptMap["-x2"] = 29;
    sptMap["+x3"] = 3; sptMap["+y3"] = 12; sptMap["-y3"] = 21; sptMap["-x3"] = 30;
    sptMap["+x4"] = 4; sptMap["+y4"] = 13; sptMap["-y4"] = 22; sptMap["-x4"] = 31;
    sptMap["+x5"] = 5; sptMap["+y5"] = 14; sptMap["-y5"] = 23; sptMap["-x5"] = 32;
    sptMap["+x6"] = 6; sptMap["+y6"] = 15; sptMap["-y6"] = 24; sptMap["-x6"] = 33;
    sptMap["+x7"] = 7; sptMap["+y7"] = 16; sptMap["-y7"] = 25; sptMap["-x7"] = 34;
    sptMap["+x8"] = 8; sptMap["+y8"] = 17; sptMap["-y8"] = 26; sptMap["-x8"] = 35;
  }
  std::map<std::string,int>::const_iterator itr = sptMap.find(spt);
  if ( itr == sptMap.end() ) return kFALSE;
  id = itr->second;
  return kTRUE;
}

Bool_t ConfigReader::getTfeId(DOMElement& elem, Int_t& id) {
  id = ChannelKey::UNDEF;
  std::string ID("ID");
  std::string tfe = Dom::getAttribute(&elem,ID);
  std::string BCAST("BCAST");
  std::string bcast("bcast");
  if ( tfe == BCAST || tfe == bcast) { 
    id = ChannelKey::BCAST;
    return kTRUE;    
  }
  static std::map<std::string,int> tfeMap;
  if ( tfeMap.size() == 0) {
    tfeMap["+x0"] = 0; tfeMap["+y0"] = 1; tfeMap["-y0"] = 2; tfeMap["-x0"] = 3;
    tfeMap["+x1"] = 4; tfeMap["+y1"] = 5; tfeMap["-y1"] = 6; tfeMap["-x1"] = 7;
    tfeMap["+x2"] = 8; tfeMap["+y2"] = 9; tfeMap["-y2"] = 10; tfeMap["-x2"] = 11;
    tfeMap["+x3"] = 12; tfeMap["+y3"] = 13; tfeMap["-y3"] = 14; tfeMap["-x3"] = 15;
    tfeMap["+x4"] = 16; tfeMap["+y4"] = 17; tfeMap["-y4"] = 18; tfeMap["-x4"] = 19;
    tfeMap["+x5"] = 20; tfeMap["+y5"] = 21; tfeMap["-y5"] = 22; tfeMap["-x5"] = 23;
  }
  std::map<std::string,int>::const_iterator itr = tfeMap.find(tfe);
  if ( itr == tfeMap.end() ) return kFALSE;
  id = itr->second;
  return kTRUE;
}


void ConfigReader::loadAcdIds() {

  ConfigBranchImpl<UShort_t>* reg = static_cast<ConfigBranchImpl<UShort_t>*>(m_config->branch("acd_id"));
  UShort_t val(0);
  UInt_t tile(0); UInt_t pmt(0);
  for ( Int_t garc(0); garc < 12; garc++ ) {
    for ( Int_t gafe(0); gafe < 18; gafe++ ) {
      AcdMap::convertToTilePmt(garc,gafe,tile,pmt);
      val = AcdMap::gemIndexFromTile(tile);
      if ( pmt != 0 ) val += 128;
      ChannelKey key(garc,gafe);  
      reg->setVal(key,val);
    }
  }
}
