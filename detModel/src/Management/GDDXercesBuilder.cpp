#include "xml/Substitute.h"
#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_Node.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_NamedNodeMap.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include "vector.h"

#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDXercesBuilder.h"


#include "detModel/Sections/GDDsection.h"
#include "detModel/Sections/GDDbox.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDaxisPos.h"
#include "detModel/Sections/GDDaxisMPos.h"
#include "detModel/Sections/GDDchoice.h"
#include "detModel/Sections/GDDidField.h"
#include "detModel/Sections/GDDseg.h"

GDDXercesBuilder::GDDXercesBuilder(char* nameFile)
{
  unsigned int iSec;

  xml::XmlParser* parser = new xml::XmlParser();

  domfile = parser->parse(nameFile);

  xml::Substitute* sub = new xml::Substitute(domfile);
  DOM_Element docElt = domfile.getDocumentElement();

  DOM_NodeList sections = docElt.getElementsByTagName(DOMString("section"));

  for (iSec = 0; iSec < sections.getLength(); iSec++) {
    DOM_Node  secNode = sections.item(iSec);
    DOM_Element& secElt = static_cast<DOM_Element &> (secNode);
    sub->execute(secElt);
  }

  GDDmanager* man = GDDmanager::getPointer();
  currentGDD = man->getGDD();
  currentGDD->setCVSid(xml::Dom::getAttribute(docElt, "CVSid"));
  currentGDD->setDTDversion(xml::Dom::getAttribute(docElt, "DTDversion"));
  
}



void GDDXercesBuilder::buildSections()
{
  unsigned int i;
  DOM_Element docElt = domfile.getDocumentElement();
  DOM_NodeList childs = docElt.getChildNodes();

  for(i=0;i<childs.getLength();i++)
    {
      std::string str = std::string(xml::Dom::transToChar(childs.item(i).getNodeName()));
      if(str == "section"){
	currentGDD->getSections()->push_back( buildSection( &(childs.item(i) )  ));
      }
    }

  currentGDD->buildVolumeMap();

  currentGDD->ResolveReferences();

}

///////////////////////////////////////////////////////////////////////////////////

GDDsection* GDDXercesBuilder::buildSection(DOM_Node* e)
{
  unsigned int i;

  DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));
  GDDsection* s = new GDDsection(
				 xml::Dom::getAttribute(el, "name"),
				 xml::Dom::getAttribute(el, "version"),
				 xml::Dom::getAttribute(el, "date"),
				 xml::Dom::getAttribute(el, "author"),
				 xml::Dom::getAttribute(el, "topVolume")
				 );
  DOM_NodeList childs = e->getChildNodes();
  for(i=0;i<childs.getLength();i++)
    {
      if (childs.item(i).getNodeType()!=Comment){
        std::string str = std::string(xml::Dom::transToChar(childs.item(i).getNodeName()));
        if(str == "box")
  	  s->addVolume(buildBox(&(childs.item(i))));
	else if(str == "choice")
	  s->addChoice(buildChoice(&(childs.item(i)))); 
        else if(str == "composition")
          s->addVolume(buildComposition( &(childs.item(i)) ));
        else if(str=="stackX"||str=="stackY"||str=="stackZ")
          s->addVolume(buildStack( &(childs.item(i)) ));
      }
    }
  return s;
}

///////////////////////////////////////////////////////////////////////////
GDDchoice* GDDXercesBuilder::buildChoice(DOM_Node* e)
{
  DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));
  GDDchoice* b = new GDDchoice(xml::Dom::getAttribute(el, "name"));

  DOM_NodeList childs = e->getChildNodes();

  for(unsigned int i=0;i<childs.getLength();i++){
      if (childs.item(i).getNodeType()!=Comment)
	{
	  DOM_NamedNodeMap attributelist=childs.item(i).getAttributes();
	  
	  std::string Mode = 
	    std::string(
			xml::Dom::transToChar(
			attributelist.getNamedItem(DOMString("mode")).getNodeValue()));

	  std::string Volume = 
	    std::string(
			xml::Dom::transToChar(
			attributelist.getNamedItem(DOMString("volume")).getNodeValue()));
	  b->addCase(Mode, Volume);
	}
  }
  return b;
}

///////////////////////////////////////////////////////////////////////////
/// \doit Add sensitive and parameters
GDDbox* GDDXercesBuilder::buildBox(DOM_Node* e)
{

  DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));
  GDDbox* b = new GDDbox(xml::Dom::getAttribute(el, "name"));

  b->setUnitLength(xml::Dom::transToChar(el.getAttribute("unitLength")));
  b->setUnitAngle(xml::Dom::transToChar(el.getAttribute("unitAngle")));
  if(el.getAttribute("X") != "0")
    b->setX(atof(xml::Dom::transToChar(el.getAttribute("X"))));
  if(el.getAttribute("Y") != "0")
    b->setY(atof(xml::Dom::transToChar(el.getAttribute("Y"))));
  if(el.getAttribute("Z") != "0")
    b->setZ(atof(xml::Dom::transToChar(el.getAttribute("Z"))));
  b->setMaterial(xml::Dom::transToChar(el.getAttribute("material")));

   if (e->hasChildNodes()){
    DOM_Node child=e->getFirstChild();
    GDDseg* s = new GDDseg;

    DOM_NamedNodeMap attributelist = child.getAttributes();
      for(unsigned int k=0;k<attributelist.getLength();k++){
       std::string NameAttr=
            std::string(xml::Dom::transToChar(attributelist.item(k).getNodeName()) );
       char* ValueAttr=xml::Dom::transToChar( attributelist.item(k).getNodeValue() );
       if(NameAttr=="axis")s->setAxis(ValueAttr);
       else if(NameAttr=="reason")s->setReason(ValueAttr);
       else if(NameAttr=="nSeg" && "1")s->setnSeg(atoi(ValueAttr));
      }//end for
     b->setSeg(s); 
   }
  return b;
}

////////////////////////////////////////////////////////////////////////

GDDcomposition* GDDXercesBuilder::buildComposition(DOM_Node* e){
  unsigned int i;

  GDDcomposition* b = new GDDcomposition();
  DOM_NamedNodeMap attributelist=e->getAttributes();

  for (i=0;i<attributelist.getLength();i++){
    if (attributelist.item(i).getNodeType()!=Comment){
      DOM_Node currentAttribute=attributelist.item(i);
      std::string value = std::string(xml::Dom::transToChar( currentAttribute.getNodeValue() ));
      std::string attributeName = 
	std::string(xml::Dom::transToChar( currentAttribute.getNodeName() ));
      
      if( attributeName=="name") b->setName(value);
      else if (attributeName=="envelope"){
	if ( value !="" )
	  b->setEnvelopeRef(value);
      }
      else if (attributeName== "parameters"){
	if ( value !="" )
	  b->setParameters(value);
      }
    }
  };

  DOM_NodeList childs = e->getChildNodes();

  for(unsigned int i=0;i<childs.getLength();i++){
      if (childs.item(i).getNodeType()!=Comment)
	b->addPosition(buildPosition( &(childs.item(i)) ));
      
  };

  return b;
}

////////////////////////////////////////////////////////////////////////////////
GDDanyPosition* GDDXercesBuilder::buildPosition(DOM_Node* e){

  std::string posType = std::string(xml::Dom::transToChar(e->getNodeName()));
  
  if (posType=="posXYZ"){
    
    GDDposXYZ* pos=new GDDposXYZ(posXYZ);
    setAttributePosition(pos, e);
       
    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));

    if(el.getAttribute("X") != "0")
      pos->setX(atof(xml::Dom::transToChar(el.getAttribute("X"))));
    if(el.getAttribute("Y") != "0")
      pos->setY(atof(xml::Dom::transToChar(el.getAttribute("Y"))));
    if(el.getAttribute("Z") != "0")
      pos->setZ(atof(xml::Dom::transToChar(el.getAttribute("Z"))));
    return pos;
  }
  else if (posType=="posRPhiZ"){
    /// \todo Add posRPhiZ
    return 0;
  }
}



/////////////////////////////////////////////////////////////////////////
GDDanyRelativePosition* GDDXercesBuilder::buildRelativePosition(DOM_Node* e){

  std::string posType = std::string(xml::Dom::transToChar(e->getNodeName()));
  //AxisPos,AxisMpos
  if (posType=="axisPos"){
      GDDaxisPos* pos=new GDDaxisPos();
      
      setAttributeRelativePosition(pos, e);
       
      DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));

      if(el.getAttribute("shift") != "0")
	pos->setShift(atof(xml::Dom::transToChar(el.getAttribute("shift"))));
      if(el.getAttribute("gap") != "0")
	pos->setGap(atof(xml::Dom::transToChar(el.getAttribute("gap"))));
    
      return pos;
      
  }
  
  else{
    GDDaxisMPos* pos=new GDDaxisMPos();

    setAttributeRelativePosition(pos, e);
       
    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));
    
    if(el.getAttribute("shift") != "0")
      pos->setShift(atof(xml::Dom::transToChar(el.getAttribute("shift"))));
    if(el.getAttribute("gap") != "0")
      pos->setGap(atof(xml::Dom::transToChar(el.getAttribute("gap"))));
    if(el.getAttribute("shift0") != "0")
      pos->setShift0(atof(xml::Dom::transToChar(el.getAttribute("shift0"))));
    if(el.getAttribute("gap0") != "0")
      pos->setGap0(atof(xml::Dom::transToChar(el.getAttribute("gap0"))));
    if(el.getAttribute("ncopy") != "0")
      pos->setNcopy(atof(xml::Dom::transToChar(el.getAttribute("ncopy"))));
    
    return pos;
    
  }
}
//////////////////////////////////////////////////////////////////////////
GDDstack* GDDXercesBuilder::buildStack(DOM_Node* e){
  stacktype st;

  std::string Stack = std::string(xml::Dom::transToChar( e->getNodeName() ));
  if(Stack=="stackX")st=sx;
  else if(Stack=="stackY")st=sy;
  else st=sz;

  GDDstack* b = new GDDstack(st);

  DOM_NamedNodeMap attributelist=e->getAttributes();
  for (unsigned int i=0;i<attributelist.getLength();i++){
    if (attributelist.item(i).getNodeType()!=Comment){
      DOM_Node currentAttribute=attributelist.item(i);
   
       std::string value = std::string(xml::Dom::transToChar( currentAttribute.getNodeValue() ));
       std::string attributeName = 
	 std::string(xml::Dom::transToChar( currentAttribute.getNodeName() ));
       
       if( attributeName=="name") b->setName(value);
       else if (attributeName=="origin"){
	 if (value=="atStart")b->setOrigin(atStart);
	 else b->setOrigin(atCentre);
       }
       else if (attributeName== "parameters")b->setParameters(value);
    }
  }

  DOM_NodeList childs = e->getChildNodes();


  for(unsigned int i=0;i<childs.getLength();i++){
    if (childs.item(i).getNodeType()!=Comment){
      GDDanyRelativePosition* p = buildRelativePosition( &(childs.item(i)));
      p->setPosDir(st);
      b->addPosition(p);
    }
  };
  
  return b;
}

void GDDXercesBuilder::setAttributePosition(GDDanyPosition* pos, DOM_Node* e)
{
  DOM_NamedNodeMap attributelistPos=e->getAttributes();

  DOM_Node node;

  node = attributelistPos.getNamedItem(DOMString("volume"));
  pos->setVolumeRef(xml::Dom::transToChar(node.getNodeValue()));

  node = attributelistPos.getNamedItem(DOMString("xrot"));
  pos->setXrot(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("yrot"));
  pos->setYrot(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("zrot"));
  pos->setZrot(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("S"));
  pos->setS(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("unitLength"));
  pos->setUnitLength(xml::Dom::transToChar(node.getNodeValue()));

  node = attributelistPos.getNamedItem(DOMString("unitAngle"));
  pos->setUnitAngle(xml::Dom::transToChar(node.getNodeValue()));


  if (e->hasChildNodes()){
    DOM_NodeList child=e->getChildNodes();
    for(unsigned int k=0; k<child.getLength() ; k++)
      {
	if(child.item(k).getNodeType() != Comment) 
	  {
	    GDDidField* field = new GDDidField;

	    DOM_NamedNodeMap attributelist = child.item(k).getAttributes();
	    for(unsigned int i=0;i<attributelist.getLength();i++){
	      
	      std::string NameAttr=
		std::string(xml::Dom::transToChar(attributelist.item(k).getNodeName()) );
	      char* ValueAttr=xml::Dom::transToChar( attributelist.item(k).getNodeValue() );
	      if(NameAttr=="name")field->setName(ValueAttr);
	      else if(NameAttr=="step" && ValueAttr!="0")field->setStep(atof(ValueAttr));
	      else if(NameAttr=="value" && ValueAttr!="0")field->setValue(atof(ValueAttr));
	    }//end for
	    pos->addIdField(field);
	  }//end if
      }
  }
}

void GDDXercesBuilder::setAttributeRelativePosition(GDDanyRelativePosition* pos, DOM_Node* e)
{
  DOM_NamedNodeMap attributelistPos=e->getAttributes();

  DOM_Node node;

  node = attributelistPos.getNamedItem(DOMString("volume"));
  pos->setVolumeRef(xml::Dom::transToChar(node.getNodeValue()));

  node = attributelistPos.getNamedItem(DOMString("dX"));
  pos->setDx(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("dY"));
  pos->setDy(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("dZ"));
  pos->setDz(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("rotation"));
  pos->setRotation(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("unitLength"));
  pos->setUnitLength(xml::Dom::transToChar(node.getNodeValue()));

  node = attributelistPos.getNamedItem(DOMString("unitAngle"));
  pos->setUnitAngle(xml::Dom::transToChar(node.getNodeValue()));


  if (e->hasChildNodes()){
    DOM_NodeList child=e->getChildNodes();
    for(unsigned int k=0; k<child.getLength() ; k++)
      {
	if(child.item(k).getNodeType() != Comment) 
	  {
	    GDDidField* field = new GDDidField;

	    DOM_NamedNodeMap attributelist = child.item(k).getAttributes();
	    for(unsigned int i=0;i<attributelist.getLength();i++){
	      
	      std::string NameAttr=
		std::string(xml::Dom::transToChar(attributelist.item(k).getNodeName()) );
	      char* ValueAttr=xml::Dom::transToChar( attributelist.item(k).getNodeValue() );
	      if(NameAttr=="name")field->setName(ValueAttr);
	      else if(NameAttr=="step" && ValueAttr!="0")field->setStep(atof(ValueAttr));
	      else if(NameAttr=="value" && ValueAttr!="0")field->setValue(atof(ValueAttr));
	    }//end for
	    pos->addIdField(field);
	  }//end if
      }
  }
 
}




