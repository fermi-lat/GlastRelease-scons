#include "xmlBase/Dom.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMCharacterData.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMTreeWalker.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#include "xmlBase/XmlParser.h"
#include "facilities/Util.h"

#include <xercesc/dom/DOMElement.hpp>

#include <string>
#include <list>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "configData/gem/TrgConfigParser.h"

TrgConfigParser::TrgConfigParser()
  :m_topElt(0){
  reset();
}

TrgConfigParser::TrgConfigParser(const char* filename) {
  XmlParser parser(true);
  reset();
  DOMDocument* doc;
  try {
    doc =parser.parse( filename );
  }
  catch (ParseException ex) {
    std::cerr << "caught exception with message " << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    assert(false);
  }
  // top level element
  m_topElt = doc->getDocumentElement();
}

TrgConfigParser::TrgConfigParser(DOMElement* el)
  :m_topElt(el){
  reset();
};

void TrgConfigParser::reset(){
  const char *tagnames[]={"WIN","AOC","SCH","TAM","TIE","ROI","periodic_rate","configuration","periodic_mode"};
  std::list<std::string> taglist;
  for (int i=0;i<9;i++)m_missingTags.push_back(tagnames[i]);
}

std::list<std::string> TrgConfigParser::remainingTags(){
  return m_missingTags;
}

bool TrgConfigParser::configComplete(){
  return (m_missingTags.size()==0);
}
  
int TrgConfigParser::parse(TrgConfig* tcf, const char* filename){
   XmlParser parser(true);
  DOMDocument* doc;
  try {
    doc =parser.parse( filename );
  }
  catch (ParseException ex) {
    std::cerr << "TrgConfigParser error: GEM xml configuration file was either not found or it could not be parsed." << std::endl;
    assert(false);
  }
  // top level element
  m_topElt = doc->getDocumentElement();
  return parse(tcf);
}
int TrgConfigParser::parse(TrgConfig* tcf, DOMElement* el){
  m_topElt=el;
  return parse(tcf);
}
  
int TrgConfigParser::parse(TrgConfig* tcf){
  assert(m_topElt);
  if(removeWhitespace(Dom::getTagName(m_topElt))!="GEM"){
    //    std::cout<<"Searching for GEM..."<<std::endl;
    std::vector<DOMElement*> defOutList;
    Dom::getChildrenByTagName( m_topElt, "GEM", defOutList );
    if (defOutList.empty()){
      std::cout<<"No GEM tag found. Exiting"<<std::endl;
      assert(0);
    }
    if (defOutList.size()>1){
      std::cout<<"More than one GEM tag found. Exiting"<<std::endl;
      assert(0);
    }  else{
      //std::cout<<"found GEM"<<std::endl;
    }
    
    m_topElt=defOutList[0];
  }
  std::vector<DOMElement*> children;
  std::string nm,ous;

  Dom::getChildrenByTagName(m_topElt, "*", children);
  for (unsigned j=0; j < children.size(); j++) {
    DOMElement* currentElement = children[j];
    if (currentElement){
      nm=removeWhitespace(Dom::getTagName(currentElement));
      
      if (nm=="periodic_rate"){
	//	ous=removeWhitespace(Dom::getTextContent(currentElement));
		//unsigned long per_rate=strtoul(ous.c_str(),0,0);
	unsigned long per_rate=content(m_topElt,"periodic_rate");
	if(per_rate&0x80000000)tcf->_tpt.setOnePPS(true);
	else tcf->_tpt.setOnePPS(false);
	tcf->_tpt.setPrescale(per_rate&0xffffff);
	m_missingTags.remove("periodic_rate");
      }
      if (nm=="configuration"){
	//ous=removeWhitespace(Dom::getTextContent(currentElement));
	//unsigned long config=strtoul(ous.c_str(),0,0);
	unsigned long config=content(m_topElt,"configuration");
	tcf->_configuration.setConfiguration(config);
	m_missingTags.remove("configuration");
      }
      if (nm=="periodic_mode"){
	unsigned long mode=content(m_topElt,"periodic_mode");
	tcf->_tpt.setFreeRunning(mode);
	m_missingTags.remove("periodic_mode");
      }
      if (nm=="WIN"){
	unsigned long win_width=content(currentElement,"window_width");
	tcf->_twp.setWindowWidth(win_width);
	m_missingTags.remove("WIN");
      }
      if (nm=="AOC"){
	unsigned long window_open_mask=content(currentElement,"window_open_mask");
	tcf->_twp.setWindowMask(window_open_mask);
	unsigned long periodic_limit=content(currentElement,"periodic_limit");
	tcf->_tpt.setLimit(periodic_limit);
	m_missingTags.remove("AOC");
      }
      if (nm=="SCH"){
	char condname[128];
	for (int i=0;i<32;i++){
	  sprintf(condname,"conditions_%02X_%02X",i*8,i*8+7);
	  unsigned long condreg=content(currentElement,condname);
	  tcf->_lut.setRegister(i,condreg);
	}
	m_missingTags.remove("SCH");
      }
      if (nm=="TAM"){
	char engname[128];
	for (int i=0;i<16;i++){
	  sprintf(engname,"engine_%x",i);
	  unsigned long engine=content(currentElement,engname);
	  tcf->_tev.setEngine(i,engine);
	}
	m_missingTags.remove("TAM");
      }
      if (nm=="TIE"){
	unsigned long reg;
	reg=content(currentElement,"towers_0_3");	
	tcf->_tdv.setTowerRegister(0,reg);
	reg=content(currentElement,"towers_4_7");	
	tcf->_tdv.setTowerRegister(1,reg);
	reg=content(currentElement,"towers_8_b");	
	tcf->_tdv.setTowerRegister(2,reg);
	reg=content(currentElement,"towers_c_f");	
	tcf->_tdv.setTowerRegister(3,reg);

	reg=content(currentElement,"tiles_000_013");	
	tcf->_tdv.setTileRegister(0,reg);
	reg=content(currentElement,"tiles_014_032");	
	tcf->_tdv.setTileRegister(1,reg);
	reg=content(currentElement,"tiles_033_NA3");	
	tcf->_tdv.setTileRegister(2,reg);
	reg=content(currentElement,"tiles_100_113");	
	tcf->_tdv.setTileRegister(3,reg);
	reg=content(currentElement,"tiles_114_NA5");	
	tcf->_tdv.setTileRegister(4,reg);
	reg=content(currentElement,"tiles_200_213");	
	tcf->_tdv.setTileRegister(5,reg);
	reg=content(currentElement,"tiles_214_NA7");	
	tcf->_tdv.setTileRegister(6,reg);
	reg=content(currentElement,"tiles_300_313");	
	tcf->_tdv.setTileRegister(7,reg);
	reg=content(currentElement,"tiles_314_NA9");	
	tcf->_tdv.setTileRegister(8,reg);
	reg=content(currentElement,"tiles_400_413");	
	tcf->_tdv.setTileRegister(9,reg);
	reg=content(currentElement,"tiles_414_NA1");	
	tcf->_tdv.setTileRegister(10,reg);
	reg=content(currentElement,"tiles_500_NA10");	
	tcf->_tdv.setTileRegister(11,reg);
	  
	reg=content(currentElement,"acd_cno");	
	tcf->_tdv.setCnoRegister(reg);
	
	reg=content(currentElement,"tower_busy");	
	tcf->_tdv.setBusyRegister(reg);
	// different versions of the xml file have different names for external 
	std::vector<DOMElement*> rList;
	Dom::getChildrenByTagName( currentElement, "external", rList );
	if(rList.size()==1) reg=content(currentElement,"external");	
	else reg=content(currentElement,"external_trg");	
	tcf->_tdv.setExternal(reg);
	m_missingTags.remove("TIE");
      }
      if (nm=="ROI"){
	unsigned long reg;
	const char* regname[]={"r_000_001", "r_002_003", "r_004_010", "r_011_012", "r_013_014", "r_020_021", "r_022_023", "r_024_030", 
		    "r_031_032", "r_033_034", "r_040_041", "r_042_043", "r_044", "r_100", "r_101_102", "r_103_104", 	 
		    "r_110_111", "r_112_113", "r_114_120", "r_121_122", "r_123_124", "r_130", "r_200", "r_201_202", 	 
		    "r_203_204", "r_210_211", "r_212_213", "r_214_220", "r_221_222", "r_223_224", "r_230", "r_300", 	 
		    "r_301_302", "r_303_304", "r_310_311", "r_312_313", "r_314_320", "r_321_322", "r_323_324", "r_330", 	 
		    "r_400", "r_401_402", "r_403_404", "r_410_411", "r_412_413", "r_414_420", "r_421_422", "r_423_424", 	 
		    "r_430", "r_500", "r_501_502", "r_503_600", "r_601_602", "r_603"};                                       
	for (int i=0;i<54;i++){
	  std::vector<DOMElement*> rList;
	  Dom::getChildrenByTagName( currentElement, regname[i], rList );	  
	  if(rList.empty()){
	    std::cout<<"Non-existant tag "<<regname[i]<<std::endl;
	    assert (0);
	  }
	  std::vector<DOMElement*> sList;
	  Dom::getChildrenByTagName( rList[0], "*", sList );	
	  if(sList.size()<2){
	    reg=content(currentElement,regname[i]);	
	  }else if(sList.size()==2){
	    std::string tag1=maketag(regname[i],1);
	    std::string tag2=maketag(regname[i],2);
	    // bug in xml file
	    if (std::string(regname[i])=="r_312_313" && Dom::getTagName(sList[1])=="t_314")tag2="t_314";
	    unsigned long c1=content(rList[0],tag1.c_str());
	    unsigned long c2=content(rList[0],tag2.c_str());
	    reg=c1|(c2<<16);
	  }else{
	    std::cout<<"Bad definition of ROI in xml file"<<std::endl;
	    assert (0);
	  }
	  tcf->_roi.setRoiRegister(i,reg);
	}
	m_missingTags.remove("ROI");
      }
    }
  }
  return 0;
}
 
const std::string TrgConfigParser::removeWhitespace(const std::string inst) const{
  std::string st(inst);
  std::string::iterator end=std::remove(st.begin(), st.end(), ' ');
  st.erase(end,st.end());
  return st;
}
unsigned long TrgConfigParser::content(DOMElement* el, const char* tag){
  std::vector<DOMElement*> subList;
  std::string ous;
  unsigned res;
  Dom::getChildrenByTagName( el, tag, subList );
  if (subList.empty() ){
    return 0;
  } else if(subList.size()>1){
    std::cout<<"More than 1 "<<tag<<" tag found. Exiting"<<std::endl;
    assert(0);
  }
  if (removeWhitespace(Dom::getTagName(subList[0]))==tag){
    std::vector<DOMElement*> rList;
    Dom::getChildrenByTagName( subList[0], "register", rList );
    if(rList.empty()){
      ous=removeWhitespace(Dom::getTextContent(subList[0]));
    }else{
      ous=removeWhitespace(Dom::getTextContent(rList[0]));
    }
    res=strtoul(ous.c_str(),0,0);
  }else {
    res=0;
  }
  return res;
}
const std::string TrgConfigParser::maketag(const char* fulltag,const int index)const{
  std::string retstr;
  if (index==1)retstr=std::string("t_")+std::string(fulltag).substr(2,3);
  else if (index==2)retstr=std::string("t_")+std::string(fulltag).substr(6,3);
  else {
    std::cout<<"maketag index has to be 1 or 2"<<std::endl;
    assert(0);
  }
  return retstr;
}
    
