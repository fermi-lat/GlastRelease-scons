// $Header$
/// Module provides methods for clients to get strip services.

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_TreeWalker.hpp"

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <memory>

#include "calibUtil/ClientObject.h"
#include "calibUtil/GenericSrv.h"

namespace calibUtil {
  

  // Initialize the data structures by parsing the XML file
  StripSrv::StripSrv(std::string xmlFileName) : m_badType(UNKNOWN_BADTYPE),
                                                m_genSrv(0)            
  {
    xml::XmlParser* parser = new xml::XmlParser();
    DOM_Document doc = parser->parse(xmlFileName.c_str());

    if (doc == DOM_Document()) {
      std::cerr << "Error parsing document" << xmlFileName << std::endl;
      return;
    }

    DOM_Element docElt = doc.getDocumentElement();
    m_genSrv = new GenericSrv(docElt);

    std::string bString = xml::Dom::getAttribute(docElt,"badType");
    if (!bString.compare("hot")) m_badType = HOT;
    else if (!bString.compare("dead")) m_badType = DEAD;
    else m_badType = UNKNOWN_BADTYPE;

    DOM_Element towerElt = xml::Dom::findFirstChildByName(docElt, "tower");

    while (towerElt != DOM_Element() ) {
      Tower tower;
      tower.id.row = atoi((xml::Dom::getAttribute(towerElt,"row")).c_str());
      tower.id.col = atoi((xml::Dom::getAttribute(towerElt,"col")).c_str());

      DOM_Element trayElt = xml::Dom::getFirstChildElement(towerElt);

      while (trayElt != DOM_Element()) {
        Tray tray;
        tray.top = tray.bot = 0;

        tray.id = atoi((xml::Dom::getAttribute(trayElt, "trayNum")).c_str());
        
        Unilayer *uni = new Unilayer;
        DOM_Element uniElt = xml::Dom::getFirstChildElement(trayElt);
        std::string which = xml::Dom::getAttribute(uniElt, "which");

        if (which == std::string("top")) tray.top = uni;
        else                             tray.bot = uni;
        fillUni(uniElt, uni);

        // Might or might not have additional unilayer child
        DOM_Element nextUniElt = xml::Dom::getSiblingElement(uniElt);
        
        if (nextUniElt != DOM_Element()) {
          Unilayer *nextUni = new Unilayer;

          if (tray.top) tray.bot = nextUni;
          else tray.top = nextUni;
          fillUni(nextUniElt, nextUni);
        }

        tower.trays.push_back(tray);
        
        trayElt = xml::Dom::getSiblingElement(trayElt);
      }
      m_towers.push_back(tower);
      towerElt = xml::Dom::getSiblingElement(towerElt);
    }
  }

  /// destructor used to deallocate memory
  StripSrv::~StripSrv(){
    delete m_genSrv;
  }
  
  /// returns the status (Hot or Dead) of the strip
  StripSrv::eBadType  StripSrv::getBadType() const {
    return m_badType;
  }
  
  /// Lists all towers with bad strips 
  void StripSrv::getBadTowers(std::vector<StripSrv::towerRC>& towerRCs) const
  {
    std::vector<Tower>::const_iterator it = m_towers.begin();
    while(it != m_towers.end() ) {

      towerRC trc;
      trc.row = it->id.row;
      trc.col = it->id.col;
      // cout << trc.row;
      // cout << trc.col;
      towerRCs.push_back(trc);
      it++;
    }
  }
  
  /// counts very bad strips for the tower specified 
  unsigned int  StripSrv::nVeryBad(towerRC towerId) const {
    unsigned nStrip = 0;

    const Tower* tower = findTower(towerId);
    if (!tower) return 0;

    const std::vector<Tray>& trays = tower->trays;

    for (unsigned int i = 0; i < trays.size(); i++) {
      Unilayer *uni = trays[i].top;
      if (uni != NULL){
        nStrip += uni->badLists[VERYBAD].size();
      }
      uni = trays[i].bot;
      if (uni != NULL){
        nStrip += uni->badLists[VERYBAD].size();
      }
    }
    return nStrip;
  }
  
  /// counts bad strips (including very bad) for the tower specified 
  unsigned int  StripSrv::nBad(towerRC towerId) const {
    unsigned nStrip = 0;
    const Tower* tower = findTower(towerId);

    if (!tower) return 0;

    const std::vector<Tray>& trays = tower->trays;

    for (unsigned int i = 0; i < trays.size(); i++) {
      Unilayer *uni = trays[i].top;
      if (uni != NULL){
        nStrip += uni->badLists[BAD].size();
        nStrip += uni->badLists[VERYBAD].size();
      }
      uni = trays[i].bot;
      if (uni != NULL){
        nStrip += uni->badLists[BAD].size();
        nStrip += uni->badLists[VERYBAD].size();
      }
    }
    return nStrip;
  }
  
  /// counts very bad strips for the tower and tray specified 
  unsigned int StripSrv::nVeryBad(towerRC towerId, unsigned int trayNum) const{

    const Tower* tower = findTower(towerId);
    if (!tower) return 0;

    const Tray* tray = findTray(*tower, trayNum);
    if (!tray) return 0;

    unsigned nStrip = 0;
    
    Unilayer *uni = tray->top;
    if (uni != NULL){
      nStrip += uni->badLists[VERYBAD].size();
    }
    uni = tray->bot;
    if (uni != NULL){
      nStrip += uni->badLists[VERYBAD].size();
    }
    return nStrip;
  }
  
  /// counts bad strips (including very bad) for the tower and tray specified 
  unsigned int  StripSrv::nBad(towerRC towerId, unsigned int trayNum) const {

    const Tower* tower = findTower(towerId);
    if (!tower) return 0;

    const Tray* tray = findTray(*tower, trayNum);
    if (!tray) return 0;

    unsigned nStrip = 0;
    
    Unilayer *uni = tray->top;
    if (uni != NULL){
      nStrip += uni->badLists[VERYBAD].size();
      nStrip += uni->badLists[BAD].size();
    }
    uni = tray->bot;
    if (uni != NULL){
      nStrip += uni->badLists[VERYBAD].size();
      nStrip += uni->badLists[BAD].size();
    }
    return nStrip;
  }
    
  /// counts very bad strips for the tower,tray and unilayer  specified 
  unsigned int  StripSrv::nVeryBad(towerRC towerId, unsigned int trayNum, 
                                   eUnilayer uni)  const {

    const Tower* tower = findTower(towerId);
    if (!tower) return 0;

    const Tray* tray = findTray(*tower, trayNum);
    if (!tray) return 0;

    unsigned nStrip = 0;

    if (uni == TOP) {
      Unilayer *pUni = tray->top;
      if  (  (pUni != NULL)  ){
        nStrip += pUni->badLists[VERYBAD].size();
      }
    }
    else if (uni == BOT) {
      Unilayer *pUni = tray->bot;
      if (pUni != NULL) {
        nStrip += pUni->badLists[VERYBAD].size();
      }
    }
    return nStrip;
  }
  
  /// counts bad strips (including very bad) for the tower, tray 
  /// and unilayer specified 
  unsigned int  StripSrv::nBad(towerRC towerId, unsigned int trayNum, 
                               eUnilayer uni) const {

    const Tower* tower = findTower(towerId);
    if (!tower) return 0;

    const Tray* tray = findTray(*tower, trayNum);
    if (!tray) return 0;

    unsigned nStrip = 0;
    if (uni == TOP) {
      Unilayer *pUni = tray->top;
      if  (  (pUni != NULL)  ){
        nStrip += pUni->badLists[BAD].size();
        nStrip += pUni->badLists[VERYBAD].size();
      }
    }
    else if (uni == BOT) {
      Unilayer *pUni = tray->bot;
      if (pUni != NULL) {
        nStrip += pUni->badLists[BAD].size();
        nStrip += pUni->badLists[VERYBAD].size();
      }
    }
    return nStrip;
  }
  
  /// Lists all very bad strips with the tower,tray and unilayer  
  
  void StripSrv::getVeryBad(towerRC towerId, unsigned int trayNum, 
                            eUnilayer uni,
                            std::vector<unsigned int>& strips) const{
    const Tower *tower = findTower(towerId);
    if (!tower) return;

    const Tray *tray = findTray(*tower, trayNum);
    if (!tray) return;

    Unilayer *pUni = 0;
    switch(uni) {
    case TOP:
      pUni = tray->top;
      break;
    case BOT:
      pUni = tray->bot;
      break;
    default:
      return;
    }
    
    strips.insert(strips.end(), pUni->badLists[VERYBAD].begin(),
                  pUni->badLists[VERYBAD].end());
    return;
  }

  
  /// Lists  bad strips (including very bad) with the tower, tray 
  /// and unilayer specified 

  void StripSrv::getBad(towerRC towerId, unsigned int trayNum, eUnilayer uni,
                   std::vector<unsigned int>&  strips) const {
    const Tower *tower = findTower(towerId);
    if (!tower) return;

    const Tray *tray = findTray(*tower, trayNum);
    if (!tray) return;

    Unilayer *pUni = 0;
    switch(uni) {
    case TOP:
      pUni = tray->top;
      break;
    case BOT:
      pUni = tray->bot;
      break;
    default:
      return;
    }
    
    strips.insert(strips.end(), pUni->badLists[VERYBAD].begin(),
                  pUni->badLists[VERYBAD].end());
    strips.insert(strips.end(), pUni->badLists[BAD].begin(),
                  pUni->badLists[BAD].end());
    return;
  }

 /// methods giving access to generic data

  /// Get instrument name
  std::string StripSrv::getInst() const{
    return m_genSrv->getInst();
  }
  
  /// Get timestamp
  std::string StripSrv::getTimestamp() const {
    return m_genSrv->getTimestamp();
  }
  
  /// Get calibration type
  std::string StripSrv::getCalType() const {
    return m_genSrv->getCalType();
  }
  
  /// Get format Version
  std::string StripSrv::getFmtVer() const {
    return m_genSrv->getFmtVer();
  } 
  
  
  /// call back method for client to access large data
  StripSrv::eRet StripSrv::traverseInfo(ClientObject *client) const {

    /* NOTE:  could also check for empty badLists and only
       call back client if nonempty.
    */
    std::vector<Tower>::const_iterator it = m_towers.begin();
    eRet userRet;
    
    while (it != m_towers.end() ){
      for (unsigned int i = 0; i < it->trays.size(); i++){
        const Tray& tray = it->trays[i];
        if (tray.top != NULL){
          StripCol* veryBad = &(tray.top->badLists[VERYBAD]);
          if (veryBad->size() > 0) {
            userRet = client->readData(it->id, tray.id, TOP, VERYBAD, veryBad);
            if (userRet != CONT) return userRet;
          }
          StripCol* bad = &(tray.top->badLists[BAD]);
          if (bad->size() > 0) {
            userRet = client->readData(it->id, tray.id, TOP, BAD, bad);
            if (userRet != CONT) return userRet;
          }

        }
        if(tray.bot != NULL){
          StripCol* veryBad = &(tray.bot->badLists[VERYBAD]);
          if (veryBad->size() > 0) {
            userRet = client->readData(it->id, tray.id, BOT, VERYBAD, veryBad);
            if (userRet != CONT) return userRet;
          }
          StripCol* bad = &(tray.bot->badLists[BAD]);
          if (bad->size() > 0) {
            userRet = client->readData(it->id, tray.id, BOT, BAD, bad);
            if (userRet != CONT) return userRet;
          }
        }
      }    
      it++;
    }
    return DONE;
  }

  // Private utilities
  void StripSrv::fillUni(const DOM_Element& uniElt, Unilayer* uni) {
    DOM_Element badElt = xml::Dom::findFirstChildByName(uniElt, "bad");
    if (badElt != DOM_Element() ) {
      fillStrips(badElt, uni->badLists[BAD]);
    }
    DOM_Element veryBadElt = xml::Dom::findFirstChildByName(uniElt, "veryBad");
    if (veryBadElt != DOM_Element() ) {
      fillStrips(veryBadElt, uni->badLists[VERYBAD]);
    }
  }

  void StripSrv::fillStrips(const DOM_Element& badElt, StripCol& list) {
    // Handle stripList elt first, if any
    DOM_Element listElt = xml::Dom::findFirstChildByName(badElt, "stripList");
    if (listElt != DOM_Element()) {
      std::string xmlList = xml::Dom::getAttribute(listElt, "strips");
      strToNum(xmlList, list);
    }

    // Then one or more stripSpan elements
    DOM_Element spanElt = xml::Dom::findFirstChildByName(badElt, "stripSpan");
    while (spanElt != DOM_Element() ) {
      std::string firstStr = xml::Dom::getAttribute(spanElt, "first");
      int first = atoi(firstStr.c_str());
      std::string lastStr = xml::Dom::getAttribute(spanElt, "last");
      int last = atoi(lastStr.c_str());
      
      if ((first >= 0) && (last >= first)) {
        for (unsigned int i = first; i <= last; i++) {
          list.push_back(i);
        }
      }
      spanElt = xml::Dom::getSiblingElement(spanElt);
    }
  }

  void StripSrv::strToNum(std::string s, std::vector<unsigned int> &v){

    std::string::iterator it = s.begin();
    
    // Maybe add something to be sure all we've got are digits
    // and blanks??

    // Skip over leading blanks, if any
    while((it != s.end()) && (*it == ' ')) it++;

    // save contiguous digits in tempStr
    while ((it != s.end()) && (*it >= '0') && (*it <= '9'))    {
      std::string tempStr;    
      tempStr += *it;
      it++;  
        
      while((it != s.end()) && (*it >= '0') && (*it <= '9')){
        tempStr += *it;
        it++;
      }

      // No more contiguous digits; skip over blanks
      while((it != s.end()) && (*it == ' ')) it++;

      // Save the converted integer
      v.push_back(atoi(tempStr.c_str()));
    }

  }

  const StripSrv::Tower* StripSrv::findTower(const towerRC& towerId) const {
    std::vector<Tower>::const_iterator it = m_towers.begin();
    while(it != m_towers.end() ) {
      if ((it->id.row == towerId.row) && (it->id.col == towerId.col)) {
        return (it);
      }
      ++it;
    }
    return 0;
  }

  const StripSrv::Tray* StripSrv::findTray(const Tower& tower, 
                                     unsigned int trayNum) const {
    const std::vector<Tray>& trays = tower.trays;
    for (unsigned int i = 0; i < trays.size(); i++) {
      if (trays[i].id == trayNum) return &trays[i];
    }

    return 0;
  }
    

}// end of namespace calibUtil
