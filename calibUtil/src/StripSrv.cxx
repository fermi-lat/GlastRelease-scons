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
  
  void StripSrv::strToNum(std::string s,std::vector<unsigned int> &v){

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

    DOMString domstr = DOMString("tower");
    DOM_NodeList list = docElt.getElementsByTagName(domstr);
    
    for (unsigned int childCount = 0; childCount < list.getLength(); 
        childCount++) {
      
      DOM_Node child_tmp = list.item(childCount);
      if(domstr.equals(child_tmp.getNodeName()) && 
         child_tmp.getParentNode() == (DOM_Node)docElt) {
        
        DOM_Element& child = static_cast<DOM_Element&>(child_tmp);
   
        std::string row = xml::Dom::getAttribute(child,"row");
        std::string col = xml::Dom::getAttribute(child,"col");
     
        Tower tower;
        tower.id.row = atoi(row.c_str());
        tower.id.col = atoi(col.c_str()); 

        // cout << std::endl << "TOWER WITH ROW " << tower.id.row << " COLUMN " << tower.id.col << std::endl;

        unsigned long whatToShow = 1 << (DOM_Node::ELEMENT_NODE -1);
        DOM_TreeWalker tree =  doc.createTreeWalker(child_tmp, whatToShow, 0, 0);
        DOM_Node dnode = tree.nextNode();
          
        while((!dnode.DOM_Node::isNull()) && 
              (DOMString("tray").equals(dnode.getNodeName()))){
            
          std::string trayNum  = xml::Dom::getAttribute(dnode,"trayNum");

          Tray tray;
          tray.id = atoi(trayNum.c_str());
          tray.topLayer = NULL;
          tray.botLayer = NULL;

          // cout << std::endl << "trayNum is" << tray.id << std::endl;
                 
          dnode = tree.nextNode();

          while((!dnode.DOM_Node::isNull()) && 
                (DOMString("unilayer").equals(dnode.getNodeName()))){

            Unilayer *uni = new Unilayer;          
            std::string which  = xml::Dom::getAttribute(dnode,"which");
          
            // cout << std::endl << "which is" << which << std::endl;

             dnode = tree.nextNode();
                
            if(DOMString("bad").equals(dnode.getNodeName())){
              // cout << std::endl << "got bad" << std::endl;
              uni->stripType = "bad";
            }
                
            if(DOMString("veryBad").equals(dnode.getNodeName())){
              // cout << std::endl << "got verybad" << std::endl;
              uni->stripType = "veryBad";
            }
    
            dnode = tree.nextNode();
                
            // Check for stripList or stripSpan
            while((!dnode.DOM_Node::isNull()) && 
                  (DOMString("stripSpan").equals(dnode.getNodeName()) || 
                   DOMString("stripList").equals(dnode.getNodeName()))){
                
              if(DOMString("stripSpan").equals(dnode.getNodeName())){
                std::string first = xml::Dom::getAttribute(dnode,"first");
                std::string last  = xml::Dom::getAttribute(dnode,"last");
             
                for(unsigned int i = atoi(first.c_str()); 
                    i <= atoi(last.c_str()); i++){
                  uni->stripCol.push_back(i);
                }
              }

              if(DOMString("stripList").equals(dnode.getNodeName())){
                std::string strips  = xml::Dom::getAttribute(dnode,"strips");
                strToNum(strips,uni->stripCol);
              }

              dnode = tree.nextNode(); 
            } // end of while stripSpan || stripList

            if (which == std::string("top")){
              tray.topLayer = uni;
            }
            else  
              if (which == std::string("bot")){
                tray.botLayer = uni;
              }
            else
              cout << "ERROR PARSING : element unilayer"; 
                
          }   // end of while unilayer
          tower.trayCol.push_back(tray);
        }     // end of while tray
        m_towers.push_back(tower);
      }       // end of if
    }         // end of for

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
  //  const std::vector<StripSrv::towerRC> StripSrv::getBadTowers(){
  void StripSrv::getBadTowers(std::vector<StripSrv::towerRC>& towerRCs) 
    const
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
    
     std::vector<Tower>::const_iterator it = m_towers.begin();
     unsigned nStrip = 0;

     while(it != m_towers.end() ) {
       if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
         for(unsigned int i = 0; i < it->trayCol.size(); i++){

           if(it->trayCol[i].topLayer != NULL){
             if(it->trayCol[i].topLayer->stripType == string("veryBad")){
               nStrip += it->trayCol[i].topLayer->stripCol.size();
             }
           }
           if(it->trayCol[i].botLayer != NULL){
             if(it->trayCol[i].botLayer->stripType == string("veryBad")){
               nStrip += it->trayCol[i].botLayer->stripCol.size();
             }
           }
         }
       }
       it++;
     }
     return nStrip;
  }

  /// counts bad strips (including very bad) for the tower specified 
  unsigned int  StripSrv::nBad(towerRC towerId) const {
   
    std::vector<Tower>::const_iterator it = m_towers.begin();
    unsigned nStrip = 0;
    
    while(it != m_towers.end() ) {
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){
          if(it->trayCol[i].topLayer != NULL){
            if(it->trayCol[i].topLayer->stripType == std::string("bad")){
               nStrip += it->trayCol[i].topLayer->stripCol.size();
            }
          }
          if(it->trayCol[i].botLayer != NULL){
            if(it->trayCol[i].botLayer->stripType == std::string("bad")){
              nStrip += it->trayCol[i].botLayer->stripCol.size();
            }
          }
        }
      }
      it++;
    }
    return nStrip;
  }
  
  /// counts very bad strips for the tower and tray specified 
  unsigned int StripSrv::nVeryBad(towerRC towerId, unsigned int trayNum) const{
     
    std::vector<Tower>::const_iterator it = m_towers.begin();
    unsigned nStrip = 0;
    
    while(it != m_towers.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if(it->trayCol[i].topLayer != NULL){
              if(it->trayCol[i].topLayer->stripType == std::string("veryBad")){
                nStrip += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if(it->trayCol[i].botLayer != NULL){
              if(it->trayCol[i].botLayer->stripType == std::string("veryBad")){
                nStrip += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }  
        }  
      }    
      it++;
    }
    
    return nStrip;
  }
  
  /// counts bad strips (including very bad) for the tower and tray specified 
  unsigned int  StripSrv::nBad(towerRC towerId, unsigned int trayNum) const {
    
    std::vector<Tower>::const_iterator it = m_towers.begin();
    unsigned nStrip = 0;
    
    while(it != m_towers.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if(it->trayCol[i].topLayer != NULL){
              if(it->trayCol[i].topLayer->stripType == std::string("bad")){
                nStrip += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if(it->trayCol[i].botLayer != NULL){
              if(it->trayCol[i].botLayer->stripType == std::string("bad")){
                nStrip += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return nStrip;
  }
    
  /// counts very bad strips for the tower,tray and unilayer  specified 
  unsigned int  StripSrv::nVeryBad(towerRC towerId, unsigned int trayNum, 
                                   eUnilayer uni)  const {

    std::vector<Tower>::const_iterator it = m_towers.begin();
    unsigned nStrip = 0;
    
    while(it != m_towers.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uni == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("veryBad")){
                nStrip += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uni == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("veryBad")){
                nStrip += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return nStrip;

  }
  
  /// counts bad strips (including very bad) for the tower, tray 
  /// and unilayer specified 
  unsigned int  StripSrv::nBad(towerRC towerId, unsigned int trayNum, 
                               eUnilayer uni) const {

    std::vector<Tower>::const_iterator it = m_towers.begin();
    unsigned nStrip = 0;
    
    while(it != m_towers.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uni == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("bad")){
                nStrip += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uni == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("bad")){
                nStrip += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return nStrip;
  }
  
  /// Lists all very bad strips with the tower,tray and unilayer  
  
  void StripSrv::getVeryBad(towerRC towerId, unsigned int trayNum, 
                            eUnilayer uni,
                            std::vector<unsigned int>& strips) const{

    std::vector<Tower>::const_iterator it = m_towers.begin();
    
    while(it != m_towers.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uni == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("veryBad")){
                strips = it->trayCol[i].topLayer->stripCol;
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uni == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("veryBad")){
                strips = it->trayCol[i].botLayer->stripCol;
              }
            }
          }    
        }
      }    
      it++;
    }
    
  }

  
  /// Lists  bad strips (including very bad) with the tower, tray 
  /// and unilayer specified 

  void StripSrv::getBad(towerRC towerId, unsigned int trayNum, eUnilayer uni,
                   std::vector<unsigned int>&  strips) const {

    std::vector<Tower>::const_iterator it = m_towers.begin();
    
    while(it != m_towers.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uni == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("bad")){
                strips = it->trayCol[i].topLayer->stripCol;
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uni == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("bad")){
                strips = it->trayCol[i].botLayer->stripCol;
              }
            }
          }    
        }
      }    
      it++;
    }
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
  void StripSrv::traverseInfo(ClientObject *client) const {
    
    std::vector<Tower>::const_iterator it = m_towers.begin();
    
    while(it != m_towers.end() ){
      for(unsigned int i = 0; i < it->trayCol.size(); i++){           
        if(it->trayCol[i].topLayer != NULL){
          client->readData(it->id,it->trayCol[i].id, TOP,
                          it->trayCol[i].topLayer->stripCol);
        }
        if(it->trayCol[i].botLayer != NULL){
          client->readData(it->id,it->trayCol[i].id, BOT,
                          it->trayCol[i].botLayer->stripCol);
        }
      }    
      it++;
    }

  }

}// end of namespace calibUtil
