#include "xmlUtil/Substitute.h"
#include "xmlUtil/Arith.h"
#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_Node.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_NamedNodeMap.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

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

#include "detModel/Constants/GDDconstants.h"

#include "detModel/Constants/GDDintConst.h"
#include "detModel/Constants/GDDfloatConst.h"
#include "detModel/Constants/GDDdoubleConst.h"
#include "detModel/Constants/GDDstringConst.h"

GDDXercesBuilder::GDDXercesBuilder()
{
}



void GDDXercesBuilder::parseFile(char* nameFile){
  unsigned int iSec;

  xml::XmlParser* parser = new xml::XmlParser();
  
  domfile = parser->parse(nameFile);

  xmlUtil::Substitute* sub = new xmlUtil::Substitute(domfile);
  DOM_Element docElt = domfile.getDocumentElement();
  DOM_Element tmp;
  DOM_Element curConst;
  
  DOM_NodeList sections = docElt.getElementsByTagName(DOMString("section"));
  
  for (iSec = 0; iSec < sections.getLength(); iSec++) {
    DOM_Node  secNode = sections.item(iSec);
    DOM_Element& secElt = static_cast<DOM_Element &> (secNode);
    sub->execute(secElt);
  }
  
  tmp = xml::Dom::findFirstChildByName(docElt, "constants" );

  tmp = xml::Dom::findFirstChildByName(tmp, "derived" );

  if (tmp != DOM_Element())
    {
      tmp = xml::Dom::findFirstChildByName(tmp, "derCategory" );
      
      while(tmp != DOM_Element())
	{
	  curConst = xml::Dom::findFirstChildByName(tmp, "const" );
	  
	  while (curConst != DOM_Element()) {
	    xmlUtil::Arith curArith(curConst);
	    double evalValue = curArith.evaluate();
	    curArith.saveValue();
	    curConst = xml::Dom::getSiblingElement(curConst);
	  }
	  tmp = xml::Dom::getSiblingElement(tmp);
	}
    }
      

  GDDmanager* man = GDDmanager::getPointer();
  currentGDD = man->getGDD();

  currentGDD->setCVSid(xml::Dom::getAttribute(docElt, "CVSid"));
  currentGDD->setDTDversion(xml::Dom::getAttribute(docElt, "DTDversion"));

  delete parser;
  delete sub;
}  


///////////////////////////////////////////////////////////////////////////////////////////////
GDDconst* GDDXercesBuilder::buildConst(DOM_Node* e){

  std::string name;
  std::string typeOfConst;
  std::string ut;
  
  DOM_NamedNodeMap attr = e->getAttributes();

  name=std::string(xml::Dom::transToChar(attr.getNamedItem(DOMString("name")).getNodeValue() ));
  ut=std::string(xml::Dom::transToChar(attr.getNamedItem(DOMString("uType")).getNodeValue() ));

  std::string elementName=std::string(xml::Dom::transToChar(e->getNodeName()));
  if (elementName=="prim"){
    typeOfConst=std::string(xml::Dom::transToChar(attr.getNamedItem(DOMString("type")).getNodeValue()));
    if (typeOfConst=="int"){
      GDDintConst* c=new GDDintConst;
      int val=atoi(xml::Dom::transToChar(attr.getNamedItem(DOMString("value")).getNodeValue()));
      c->setName(name);
      c->setConstMeaning(ut);
      c->setValue(val);
      c->setNote(std::string(xml::Dom::transToChar(e->getFirstChild().getNodeValue())));
      return c;
    }
    else if (typeOfConst=="float"){
      GDDfloatConst* c=new GDDfloatConst;
      float val=atof(xml::Dom::transToChar(attr.getNamedItem(DOMString("value")).getNodeValue()));
      c->setName(name);
      c->setConstMeaning(ut);
      c->setValue(val);
      c->setNote(std::string(xml::Dom::transToChar(e->getFirstChild().getNodeValue())));
      return c;
    }
    else if (typeOfConst=="double"){
      GDDdoubleConst* c=new GDDdoubleConst;
      double val=atof(xml::Dom::transToChar(attr.getNamedItem(DOMString("value")).getNodeValue()));
      c->setName(name);
      c->setConstMeaning(ut);
      c->setValue(val);
      c->setNote(std::string(xml::Dom::transToChar(e->getFirstChild().getNodeValue())));
      return c;
    }
    else if (typeOfConst=="string"){
      GDDstringConst* c=new GDDstringConst;
      std::string val=std::string(xml::Dom::transToChar(attr.getNamedItem(DOMString("value")).getNodeValue()));
      c->setName(name);
      c->setConstMeaning(ut);
      c->setValue(val);
      c->setNote(std::string(xml::Dom::transToChar(e->getFirstChild().getNodeValue())));
      return c;
    }
    else return 0;
  }//end if
  else{
    //e is a const element    
    GDDdoubleConst* c=new GDDdoubleConst;
    DOM_Element n;

    double val=atof(xml::Dom::transToChar(attr.getNamedItem(DOMString("value")).getNodeValue()));

    c->setName(name);
    c->setConstMeaning(ut);
    c->setValue(val);
    
    n = static_cast<DOM_Element&>(*e);
    DOM_NodeList nodelist = n.getElementsByTagName(DOMString("notes"));

    if (nodelist.getLength())
      c->setNote(std::string(xml::Dom::transToChar(nodelist.item(0).getFirstChild().getNodeValue())));
    
    return c;
  }//end else
}

///////////////////////////////////////////////////////////////////////////////////////////////
void GDDXercesBuilder::buildConstants(){
  unsigned int i,j;
  DOM_Element docElt = domfile.getDocumentElement();
  DOM_NodeList child = docElt.getElementsByTagName(DOMString("constants"));
  if (child.getLength()){
    GDDconstants* ConstantsBranch = new GDDconstants();
    
    DOM_NodeList childs = child.item(0).getChildNodes();
    // version primary ?derived 
    for(i=0;i<childs.getLength();i++){
      if(childs.item(i).getNodeType() != Comment)
	{
	  std::string str = std::string(xml::Dom::transToChar(childs.item(i).getNodeName()));
	  if(str == "version"){
	    std::string s1,s2;
	    DOM_NamedNodeMap attr=childs.item(i).getAttributes ();
	    s1=std::string(xml::Dom::transToChar(attr.getNamedItem(DOMString("major")).getNodeValue()));
	    s2=std::string(xml::Dom::transToChar(attr.getNamedItem(DOMString("minor")).getNodeValue()));
	    ConstantsBranch->setVersion(s1,s2);
	  }
	  else if(str == "primary" || str == "derived"){
	    DOM_NodeList child = childs.item(i).getChildNodes();
	    if (child.getLength()!=0){
	      for(j=0;j<child.getLength();j++){
		GDDconstCategory* cat= new GDDconstCategory;
		std::string s1;
		DOM_NamedNodeMap attrCat=child.item(j).getAttributes();
		s1=std::string(xml::Dom::transToChar(attrCat.getNamedItem(DOMString("name")).getNodeValue()));
		cat->setName(s1);
		if(attrCat.getLength()>1){ //if the actual category has more than one attribute it is derived
		  std::string s2;
		  s2=std::string(xml::Dom::transToChar(attrCat.getNamedItem(DOMString("save")).getNodeValue()));
		  if (s2=="true")
		    cat->setSave(true);
		  cat->setPrimary(false);
		}
		DOM_Node over=child.item(j).getFirstChild().getFirstChild();
		s1=std::string(xml::Dom::transToChar( over.getNodeValue()));
		cat->setOverview(s1);
		child.item(j).removeChild(child.item(j).getFirstChild());
		//elt is the list of prim/const
		DOM_NodeList elt=child.item(j).getChildNodes();
		unsigned int k;
		for(k=0;k<elt.getLength();k++){
		  if (elt.item(k).getNodeType()!=Comment){
		    cat->addConstant(buildConst( &elt.item(k)) );
		  }//end if 
		}//end for 
		//Insert a new category in GDDconstast object
		ConstantsBranch->addConstantCategory(cat);
	      
	      }//end for prim and derived cateory 
	    }
	  }  
	}
    }//end for                                           
    currentGDD->setConstants(ConstantsBranch);
    currentGDD->buildConstantsMap();
  }//end if
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
	currentGDD->addSection(buildSection(&(childs.item(i))));
      }
    }
  currentGDD->buildVolumeMap();

  currentGDD->ResolveReferences();

  currentGDD->buildBoundingBoxes();

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

  GDDmanager* gddManager = GDDmanager::getPointer();

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
	  
	  /// This initialize the mode
	  if(gddManager->getMode() == "")
	    gddManager->setMode(Mode);
	}
  }
  return b;
}

///////////////////////////////////////////////////////////////////////////
/// \doit Add parameters
GDDbox* GDDXercesBuilder::buildBox(DOM_Node* e)
{

  DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));
  GDDbox* b = new GDDbox(xml::Dom::getAttribute(el, "name"));

  // b->setUnitLength(xml::Dom::transToChar(el.getAttribute("unitLength")));
  // b->setUnitAngle(xml::Dom::transToChar(el.getAttribute("unitAngle")));
  if(el.getAttribute("X") != "0")
    b->setX(atof(xml::Dom::transToChar(el.getAttribute("X"))));
  if(el.getAttribute("Y") != "0")
    b->setY(atof(xml::Dom::transToChar(el.getAttribute("Y"))));
  if(el.getAttribute("Z") != "0")
    b->setZ(atof(xml::Dom::transToChar(el.getAttribute("Z"))));
  b->setMaterial(xml::Dom::transToChar(el.getAttribute("material")));

  if (std::string((xml::Dom::transToChar(el.getAttribute("sensitive")))) == "true")
    b->setSensitive(1);
  else 
    b->setSensitive(0);

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
  
  for(i=0;i<childs.getLength();i++){
    if (childs.item(i).getNodeType()!=Comment)
      b->addPosition(buildPosition( &(childs.item(i)) ));
      
  };

  return b;
}

//////////////////////////////////////////////////////////////////////////
GDDstack* GDDXercesBuilder::buildStack(DOM_Node* e){
  GDDstack::axisDir st;
  unsigned int i;

  std::string Stack = std::string(xml::Dom::transToChar( e->getNodeName() ));
  if(Stack=="stackX")st=GDDstack::xDir;
  else if(Stack=="stackY")st=GDDstack::yDir;
  else st=GDDstack::zDir;

  GDDstack* b = new GDDstack(st);

  DOM_NamedNodeMap attributelist=e->getAttributes();
  for (i=0;i<attributelist.getLength();i++){
    if (attributelist.item(i).getNodeType()!=Comment){
      DOM_Node currentAttribute=attributelist.item(i);
   
       std::string value = std::string(xml::Dom::transToChar( currentAttribute.getNodeValue() ));
       std::string attributeName = 
	 std::string(xml::Dom::transToChar( currentAttribute.getNodeName() ));
       
       if( attributeName=="name") b->setName(value);
       else if (attributeName=="origin"){
	 if (value=="atStart")b->setOrigin(GDDstack::atStart);
	 else b->setOrigin(GDDstack::atCentre);
       }
       else if (attributeName== "parameters")b->setParameters(value);
    }
  }

  DOM_NodeList childs = e->getChildNodes();
  
  
  for(i=0;i<childs.getLength();i++){
    if ( (childs.item(i)).getNodeType()!=Comment){
      GDDstackedPos* p = buildRelativePosition( &(childs.item(i)));
      p->setAxisDir((GDDstackedPos::axisDir)st);
      b->addPosition(p);
    }
  };
  
  return b;
}

////////////////////////////////////////////////////////////////////////////////
GDDsinglePos* GDDXercesBuilder::buildPosition(DOM_Node* e){

  std::string posType = std::string(xml::Dom::transToChar(e->getNodeName()));
  
  if (posType=="posXYZ"){
    
    GDDposXYZ* pos=new GDDposXYZ();
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
  else  return 0; // default
}



/////////////////////////////////////////////////////////////////////////
GDDstackedPos* GDDXercesBuilder::buildRelativePosition(DOM_Node* e){

  std::string posType = std::string(xml::Dom::transToChar(e->getNodeName()));
  //AxisPos,AxisMpos
  if (posType=="axisPos"){
      GDDaxisPos* pos=new GDDaxisPos();
      
      setAttributeRelativePosition(pos, e);
       
    
      return pos;
      
  }
  
  else{
    GDDaxisMPos* pos=new GDDaxisMPos();

    setAttributeRelativePosition(pos, e);
       
    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(*e));

    if(el.getAttribute("shift0") != "0")
      pos->setShift0(atof(xml::Dom::transToChar(el.getAttribute("shift0"))));
    if(el.getAttribute("gap0") != "0")
      pos->setGap0(atof(xml::Dom::transToChar(el.getAttribute("gap0"))));
    if(el.getAttribute("ncopy") != "0")
      pos->setNcopy(atof(xml::Dom::transToChar(el.getAttribute("ncopy"))));
    
    return pos;
    
  }
}

void GDDXercesBuilder::setAttributePosition(GDDsinglePos* pos, DOM_Node* e)
{
  DOM_NamedNodeMap attributelistPos=e->getAttributes();

  DOM_Node node;

  node = attributelistPos.getNamedItem(DOMString("volume"));
  pos->setVolumeRef(xml::Dom::transToChar(node.getNodeValue()));

  node = attributelistPos.getNamedItem(DOMString("xrot"));
  pos->setXRot(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("yrot"));
  pos->setYRot(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("zrot"));
  pos->setZRot(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("S"));
  pos->setS(atof(xml::Dom::transToChar(node.getNodeValue())));

  //  node = attributelistPos.getNamedItem(DOMString("unitLength"));
  // pos->setUnitLength(xml::Dom::transToChar(node.getNodeValue()));

  // node = attributelistPos.getNamedItem(DOMString("unitAngle"));
  // pos->setUnitAngle(xml::Dom::transToChar(node.getNodeValue()));


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
		std::string(xml::Dom::transToChar(attributelist.item(i).getNodeName()) );
	      char* ValueAttr=xml::Dom::transToChar( attributelist.item(i).getNodeValue() );
	      if(NameAttr=="name")field->setName(ValueAttr);
	      else if(NameAttr=="step")field->setStep(atof(ValueAttr));
	      else if(NameAttr=="value")field->setValue(atof(ValueAttr));
	    }//end for
	    pos->addIdField(field);
	  }//end if
      }
  }
}

void GDDXercesBuilder::setAttributeRelativePosition(GDDstackedPos* pos, DOM_Node* e)
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

  // node = attributelistPos.getNamedItem(DOMString("unitLength"));
  // pos->setUnitLength(xml::Dom::transToChar(node.getNodeValue()));

  // node = attributelistPos.getNamedItem(DOMString("unitAngle"));
  // pos->setUnitAngle(xml::Dom::transToChar(node.getNodeValue()));

  node = attributelistPos.getNamedItem(DOMString("gap"));
  pos->setGap(atof(xml::Dom::transToChar(node.getNodeValue())));

  node = attributelistPos.getNamedItem(DOMString("shift"));
  pos->setShift(atof(xml::Dom::transToChar(node.getNodeValue())));

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
		std::string(xml::Dom::transToChar(attributelist.item(i).getNodeName()) );
	      char* ValueAttr=xml::Dom::transToChar( attributelist.item(i).getNodeValue() );
	      if(NameAttr=="name")field->setName(ValueAttr);
	      else if(NameAttr=="step")field->setStep(atof(ValueAttr));
	      else if(NameAttr=="value")field->setValue(atof(ValueAttr));
	    }//end for
	    pos->addIdField(field);
	  }//end if
      }
  }
 
}




