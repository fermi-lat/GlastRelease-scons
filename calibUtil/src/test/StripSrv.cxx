// $Header$
/// Module provides methods for clients to get strip services.

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_TreeWalker.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>

#include "calibUtil/StripSrv.h"

namespace calibUtil {

  /// Initialize the data structures by parsing the XML file
  StripSrv::StripSrv(std::string xmlFileName){
    
    int num[100];
    int globalCount = 0;
    
    xml::XmlParser* parser = new xml::XmlParser();
    DOM_Document doc = parser->parse(xmlFileName.c_str());
    
    if (doc != 0) {  // successful
      std::cout << "Document successfully parsed" << std::endl;
    }
    else {     
      std::cout << "Error parsing document" << std::endl;
      exit(1);    
    }
    
    DOMString domstr = DOMString("tower");
    DOM_Element docElt = doc.getDocumentElement();
    DOM_NodeList list = docElt.getElementsByTagName(domstr);
    
    for(int childCount = 0; childCount < list.getLength(); childCount++) {
      
      DOM_Node child_tmp = list.item(childCount);
      if(domstr.equals(child_tmp.getNodeName()) && 
         child_tmp.getParentNode() == (DOM_Node)docElt) {
        
        DOM_Element& child = static_cast<DOM_Element&>(child_tmp);
        
        /* print the tower */
        std::string row = xml::Dom::getAttribute(child,"row");
        std::string col = xml::Dom::getAttribute(child,"col");
        cout << std::endl << "TOWER WITH ROW " << row << " COLUMN " << col << std::endl;
        
        



        /*..........................................*/
        




          /*

	  // print bad strips
          unsigned long whatToShow = 1 << (DOM_Node::ELEMENT_NODE -1);
	  DOM_TreeWalker tree =  doc.createTreeWalker(child_tmp, whatToShow, 0, 0);
	  //tree.setCurrentNode(child_tmp);

	  DOM_Node dnode = tree.nextNode();
          // DOMString domstr1 = DOMString("stripList");
         
          while(!dnode.DOM_Node::isNull()){
            // cout << "Checking for Strip List ";
	    if(DOMString("stripList").equals(dnode.getNodeName())){
	      std::string strips  = xml::Dom::getAttribute(dnode,"strips");
              *out << std::endl << "bad strips are " << strips << std::endl;

              // char* pcs = strips.c_str();
              char str[100];
 
              memcpy(str,strips.c_str(),strips.length());
              str[strips.length()] = '\0';
 
              *out << std::endl << "string is" << str << std::endl;

             // counting num of strips 
             char *tmpString = str;
             int  stripCount = 0; */

             /*  
             while(*tmpString != '\0')
             {
               while((*tmpString >= '0') && (*tmpString <= '9'))
               {
                 stripCount++;
                 tmpString++; 
                 while((*tmpString >= '0') && (*tmpString <= '9'))
                   tmpString++;
                 while(*tmpString == ' ')
                   tmpString++;
               }

             }

  *out << std::endl << "num of bad strips are " << stripCount<< std::endl;
             globalCount += stripCount;
             */

          /*
            strToNum(str,num);
            *out << std::endl << "trying function call, count is" << num[0] << std::endl;
           

            for(int i = 1; i <= num[0]; i++){
              cout << num[i] << "  "; 
            }
            
            } */  /* end of stripList */
          /*
              // cout << "Checking for Strip Span";
	    if(DOMString("stripSpan").equals(dnode.getNodeName())){
	      std::string first  = xml::Dom::getAttribute(dnode,"first");
              std::string last  = xml::Dom::getAttribute(dnode,"last");
              //              cout << " i am here";               
          
    int sCount = atoi(last.c_str()) -  atoi(first.c_str()) + 1;
              *out << std::endl << "num of bad strips in this Span are" << sCount<< std::endl;   
           
              // strToNum(str,num);
             globalCount += sCount;
            } 
 
            dnode = tree.nextNode();
          }
	  
 *out << std::endl << "TOTAL NUMBER OF BAD STRIPS ARE " << globalCount<< std::endl;	
	 
 */

 /*
	    DOM_Element stripList = xml::Dom::findFirstChildByName(child, "stripList");    if(!stripList.DOM_Node::isNull()){    
	    std::string strips  = xml::Dom::getAttribute(stripList,"strips");
	    *out << std::endl << "bad strips are " << strips << std::endl;
	    } 
	    
 */     
	  
          //	}
          // }
      
          // }



        /*.....................................................*/


      }
      
    }
    
  }
  

  /// return the status (Hot or Dead) of the strip
  std::string  StripSrv::getBadType(){
    std::string line("testing");
    return line;
  }
  
  /// List all towers with bad strips 
  std::vector<towerRC> StripSrv::getBadTowers(){
  }
  
  /// count for very bad strips for the tower specified 
  unsigned int  StripSrv::countVeryBad(towerRC towerId){

  }

  /// count for  bad strips (including very bad) for the tower specified 
  unsigned int  StripSrv::countBad(towerRC towerId){

  }
  
  /// count for very bad strips for the tower and tray specified 
  unsigned int  StripSrv::countVeryBad(towerRC towerId, unsigned int trayNum){



  }
  
  /// count for  bad strips (including very bad) for the tower and tray 
    /// specified 
   
  unsigned int  StripSrv::countBad(towerRC towerId, unsigned int trayNum){

  }
    
  /// count for very bad strips for the tower,tray and unilayer  specified 
  unsigned int  StripSrv::countVeryBad(towerRC towerId, unsigned int trayNum, std::string uniLayer){


  }
  
  /// count for  bad strips (including very bad) for the tower, tray 
  /// and unilayer specified 
  unsigned int  StripSrv::countBad(towerRC towerId, unsigned int trayNum,std::string uniLayer){


  }
  
  /// List all very bad strips with the tower,tray and unilayer  
  std::vector<unsigned int>  StripSrv::getVeryBad(towerRC towerId, unsigned int trayNum, std::string uniLayer){



  }

  
  /// List  bad strips (including very bad) with the tower, tray 
  /// and unilayer specified 
  std::vector<unsigned int>  StripSrv::getBad(towerRC towerId, unsigned int trayNum, std::string uniLayer){


  }

  
}/// end of namespace calibUtil






















