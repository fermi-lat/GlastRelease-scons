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
  
  /// convert stripList in string format to a vector of strip numbers
  StripSrv::strToNum(std::string s,std::vector<unsigned int> &v){

    std::string::iterator it = s.begin();
    
    // Accessing characters using iterator    
    while((it != s.end()) && (*it >= '0') && (*it <= '9'))
    {
        std::string tempStr;    
        tempStr += *it;
        it++;  
        
        while((it != s.end()) && (*it >= '0') && (*it <= '9')){
          tempStr += *it;
          it++;
        }
        while((it != s.end()) && (*it == ' '))

          it++;        
        v.push_back(atoi(tempStr.c_str()));
    }

  }


  // Initialize the data structures by parsing the XML file
  StripSrv::StripSrv(std::string xmlFileName){
    
    xml::XmlParser* parser = new xml::XmlParser();
    DOM_Document doc = parser->parse(xmlFileName.c_str());
    
    if (doc != 0) {  // successful
      std::cout << "Document successfully parsed" << std::endl;
    }
    else {     
      std::cout << "Error parsing document" << std::endl;
      exit(1);    
    }

    DOM_Element docElt = doc.getDocumentElement();
    genSrv = new GenericSrv(docElt);

    badType = xml::Dom::getAttribute(docElt,"badType");
    DOMString domstr = DOMString("tower");
    DOM_NodeList list = docElt.getElementsByTagName(domstr);
    
    for(unsigned int childCount = 0; childCount < list.getLength(); 
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

        cout << std::endl << "TOWER WITH ROW " << tower.id.row << " COLUMN " << tower.id.col << std::endl;

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

          cout << std::endl << "trayNum is" << tray.id << std::endl;
                 
          dnode = tree.nextNode();

          while((!dnode.DOM_Node::isNull()) && 
                (DOMString("unilayer").equals(dnode.getNodeName()))){

            UniLayer *uniL = new UniLayer;          
            std::string which  = xml::Dom::getAttribute(dnode,"which");
          
            cout << std::endl << "which is" << which << std::endl;

             dnode = tree.nextNode();
                
            if(DOMString("bad").equals(dnode.getNodeName())){
              cout << std::endl << "got bad" << std::endl;
              uniL->stripType = "bad";
            }
                
            if(DOMString("veryBad").equals(dnode.getNodeName())){
              cout << std::endl << "got verybad" << std::endl;
              uniL->stripType = "veryBad";
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
                  uniL->stripCol.push_back(i);
                }
              }

              if(DOMString("stripList").equals(dnode.getNodeName())){
                std::string strips  = xml::Dom::getAttribute(dnode,"strips");
                strToNum(strips,uniL->stripCol);
              }

              dnode = tree.nextNode(); 
            } // end of while stripSpan || stripList

            if(which == std::string("top")){
              tray.topLayer = uniL;
            }
            else  
              if(which == std::string("bot")){
                tray.botLayer = uniL;
              }
            else
              cout << "ERROR PARSING : element unilayer"; 
                
          }   // end of while unilayer
          tower.trayCol.push_back(tray);
        }     // end of while tray
        towerCol.push_back(tower);
      }       // end of if
    }         // end of for

  }

  /// destructor used to deallocate memory
  StripSrv::~StripSrv(){
    delete genSrv;
  }
  
  /// returns the status (Hot or Dead) of the strip
  std::string  StripSrv::getBadType(){
    return badType;
  }
  
  /// Lists all towers with bad strips 
  const std::vector<towerRC> StripSrv::getBadTowers(){
    
     std::vector<towerRC> tv;
     std::vector<Tower>::iterator it = towerCol.begin();

     while(it != towerCol.end() ) {

       towerRC trc;
       trc.row = it->id.row;
       trc.col = it->id.col;
       cout << trc.row;
       cout << trc.col;
       tv.push_back(trc);
       it++;
     }
     return tv;
  }
  
  /// counts very bad strips for the tower specified 
  unsigned int  StripSrv::countVeryBad(towerRC towerId){
    
     std::vector<Tower>::iterator it = towerCol.begin();
     unsigned stripCount = 0;

     while(it != towerCol.end() ) {
       if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
         for(unsigned int i = 0; i < it->trayCol.size(); i++){

           if(it->trayCol[i].topLayer != NULL){
             if(it->trayCol[i].topLayer->stripType == string("veryBad")){
               stripCount += it->trayCol[i].topLayer->stripCol.size();
             }
           }
           if(it->trayCol[i].botLayer != NULL){
             if(it->trayCol[i].botLayer->stripType == string("veryBad")){
               stripCount += it->trayCol[i].botLayer->stripCol.size();
             }
           }
         }
       }
       it++;
     }
     return stripCount;
  }

  /// counts bad strips (including very bad) for the tower specified 
  unsigned int  StripSrv::countBad(towerRC towerId){
   
    std::vector<Tower>::iterator it = towerCol.begin();
    unsigned stripCount = 0;
    
    while(it != towerCol.end() ) {
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){
          if(it->trayCol[i].topLayer != NULL){
            if(it->trayCol[i].topLayer->stripType == std::string("bad")){
               stripCount += it->trayCol[i].topLayer->stripCol.size();
            }
          }
          if(it->trayCol[i].botLayer != NULL){
            if(it->trayCol[i].botLayer->stripType == std::string("bad")){
              stripCount += it->trayCol[i].botLayer->stripCol.size();
            }
          }
        }
      }
      it++;
    }
    return stripCount;
  }
  
  /// counts very bad strips for the tower and tray specified 
  unsigned int  StripSrv::countVeryBad(towerRC towerId, unsigned int trayNum){
     
    std::vector<Tower>::iterator it = towerCol.begin();
    unsigned stripCount = 0;
    
    while(it != towerCol.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if(it->trayCol[i].topLayer != NULL){
              if(it->trayCol[i].topLayer->stripType == std::string("veryBad")){
                stripCount += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if(it->trayCol[i].botLayer != NULL){
              if(it->trayCol[i].botLayer->stripType == std::string("veryBad")){
                stripCount += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }  
        }  
      }    
      it++;
    }
    
    return stripCount;
  }
  
  /// counts bad strips (including very bad) for the tower and tray specified 
  unsigned int  StripSrv::countBad(towerRC towerId, unsigned int trayNum){
    
    std::vector<Tower>::iterator it = towerCol.begin();
    unsigned stripCount = 0;
    
    while(it != towerCol.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if(it->trayCol[i].topLayer != NULL){
              if(it->trayCol[i].topLayer->stripType == std::string("bad")){
                stripCount += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if(it->trayCol[i].botLayer != NULL){
              if(it->trayCol[i].botLayer->stripType == std::string("bad")){
                stripCount += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return stripCount;
  }
    
  /// counts very bad strips for the tower,tray and unilayer  specified 
  unsigned int  StripSrv::countVeryBad(towerRC towerId, 
                                       unsigned int trayNum, uniL uniLayer){

    std::vector<Tower>::iterator it = towerCol.begin();
    unsigned stripCount = 0;
    
    while(it != towerCol.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uniLayer == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("veryBad")){
                stripCount += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uniLayer == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("veryBad")){
                stripCount += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return stripCount;

  }
  
  /// counts bad strips (including very bad) for the tower, tray 
  /// and unilayer specified 
  unsigned int  StripSrv::countBad(towerRC towerId, unsigned int trayNum, 
                                   uniL uniLayer){

    std::vector<Tower>::iterator it = towerCol.begin();
    unsigned stripCount = 0;
    
    while(it != towerCol.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uniLayer == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("bad")){
                stripCount += it->trayCol[i].topLayer->stripCol.size();
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uniLayer == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("bad")){
                stripCount += it->trayCol[i].botLayer->stripCol.size();
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return stripCount;
  }
  
  /// Lists all very bad strips with the tower,tray and unilayer  
  std::vector<unsigned int>  StripSrv::getVeryBad(towerRC towerId, unsigned int trayNum, uniL uniLayer){

    std::vector<Tower>::iterator it = towerCol.begin();
    std::vector<unsigned int> stripList;
    
    while(it != towerCol.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uniLayer == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("veryBad")){
                stripList = it->trayCol[i].topLayer->stripCol;
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uniLayer == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("veryBad")){
                stripList = it->trayCol[i].botLayer->stripCol;
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return stripList;

  }

  
  /// Lists  bad strips (including very bad) with the tower, tray 
  /// and unilayer specified 
  std::vector<unsigned int>  StripSrv::getBad(towerRC towerId, unsigned int trayNum, uniL uniLayer){

    std::vector<Tower>::iterator it = towerCol.begin();
    std::vector<unsigned int> stripList;
    
    while(it != towerCol.end() ){
      if((it->id.row == towerId.row) && (it->id.col == towerId.col)){
        for(unsigned int i = 0; i < it->trayCol.size(); i++){ 
          if(it->trayCol[i].id == trayNum){ 
            if((it->trayCol[i].topLayer != NULL) && 
               (uniLayer == TOP)){
              if(it->trayCol[i].topLayer->stripType == std::string("bad")){
                stripList = it->trayCol[i].topLayer->stripCol;
              }
            }
            if((it->trayCol[i].botLayer != NULL) && 
               (uniLayer == BOT)){
              if(it->trayCol[i].botLayer->stripType == std::string("bad")){
                stripList = it->trayCol[i].botLayer->stripCol;
              }
            }
          }    
        }
      }    
      it++;
    }
    
    return stripList;
  }

 /// methods giving access to generic data

    /// Get instrument name
    std::string StripSrv::getInst(){
      return genSrv->getInst();
    }
    
    /// Get timestamp
    std::string StripSrv::getTimestamp(){
      return genSrv->getTimestamp();
    }

    /// Get calibration type
    std::string StripSrv::getCalType(){
      return genSrv->getCalType();
    }
    
    /// Get format Version
    std::string StripSrv::getFmtVer(){
      return genSrv->getFmtVer();
    } 


  /// call back method for client to access large data
  void StripSrv::traverseInfo(ClientObject *client){
    
    std::vector<Tower>::iterator it = towerCol.begin();
    
    while(it != towerCol.end() ){
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












