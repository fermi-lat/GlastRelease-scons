#include "xmlUtil/Substitute.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/docMan/GDDDocMan.h"

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

#include "detModel/Management/Manager.h"
#include "detModel/Management/XercesBuilder.h"
#include "detModel/Management/DMDocClient.h"
#include "detModel/Sections/Section.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Tube.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Sections/Choice.h"
#include "detModel/Sections/IdField.h"
#include "detModel/Sections/Seg.h"
#include "detModel/Constants/Constants.h"
#include "detModel/Constants/IntConst.h"
#include "detModel/Constants/FloatConst.h"
#include "detModel/Constants/DoubleConst.h"
#include "detModel/Constants/StringConst.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Materials/Element.h"
#include "detModel/Materials/Composite.h"

namespace detModel{

  XercesBuilder::XercesBuilder()
  {
    m_docClient = new DMDocClient();
    
    xmlUtil::GDDDocMan* pGDDMan = xmlUtil::GDDDocMan::getPointer();
    pGDDMan->regClient("constants", m_docClient);
    pGDDMan->regClient("section", m_docClient);    
    pGDDMan->regClient("materials", m_docClient);    
  }
  
  XercesBuilder::~XercesBuilder()
  {
    delete m_docClient;
  }

  void XercesBuilder::parseFile(std::string nameFile){
//    unsigned int iSec;
    xmlUtil::GDDDocMan* pGDDMan = xmlUtil::GDDDocMan::getPointer();
    
    pGDDMan->parse(nameFile);
    
    /// We start detModel stuff retriving the manager
    Manager* man = Manager::getPointer();

    /// We point to the current Gdd
    currentGdd = man->getGdd();
  
    /// Set some info on the parsed file
    //    currentGdd->setCVSid(xml::Dom::getAttribute(docElt, "CVSid"));
    // currentGdd->setDTDversion(xml::Dom::getAttribute(docElt, "DTDversion"));
  }  


  /*
   *  This method build a constant
   */
  Const* XercesBuilder::buildConst(DOM_Node e){
    std::string name;
    std::string typeOfConst;
    std::string ut;
  
    DOM_NamedNodeMap attr = e.getAttributes();

    name=std::string(xml::Dom::transToChar
                     (attr.getNamedItem(DOMString("name")).getNodeValue() ));
    ut=std::string(xml::Dom::transToChar
                   (attr.getNamedItem(DOMString("uType")).getNodeValue() ));

    std::string elementName=std::string(xml::Dom::transToChar
                                        (e.getNodeName()));
    if (elementName=="prim"){
      typeOfConst=
        std::string(xml::Dom::transToChar
                    (attr.getNamedItem(DOMString("type")).getNodeValue()));
      Const* c;
      if (typeOfConst=="int"){
	c=new IntConst;
	int val=atoi(xml::Dom::transToChar
                     (attr.getNamedItem(DOMString("value")).getNodeValue()));
	((IntConst*)c)->setValue(val);
      }
      else if (typeOfConst=="float"){
	c=new FloatConst;
	float val=atof(xml::Dom::transToChar
                       (attr.getNamedItem(DOMString("value")).getNodeValue()));
	((FloatConst*)c)->setValue(val);
      }
      else if (typeOfConst=="double"){
	c=new DoubleConst;
	double val=
          atof(xml::Dom::transToChar
               (attr.getNamedItem(DOMString("value")).getNodeValue()));
	((DoubleConst*)c)->setValue(val);
      }
      else if (typeOfConst=="string"){
	c=new StringConst;
	std::string val=
          std::string(xml::Dom::transToChar
                      (attr.getNamedItem(DOMString("value")).getNodeValue()));
	((StringConst*)c)->setValue(val);
      }
      else return 0;
      c->setName(name);
      c->setConstMeaning(ut);
      c->setNote(std::string(xml::Dom::transToChar
                             (e.getFirstChild().getNodeValue())));
      return c;
    }//end if
    else{
      //e is a const element    
      DoubleConst* c=new DoubleConst;
      DOM_Element n;

      double val=atof(xml::Dom::transToChar
                      (attr.getNamedItem(DOMString("value")).getNodeValue()));

      c->setName(name);
      c->setConstMeaning(ut);
      c->setValue(val);
    
      n = static_cast<DOM_Element&>(e);
      DOM_NodeList nodelist = n.getElementsByTagName(DOMString("notes"));

      if (nodelist.getLength()) {
	c->setNote(std::string
                   (xml::Dom::transToChar
                    (nodelist.item(0).getFirstChild().getNodeValue())));
      }
      return c;
    }//end else
  }

  /// This methods build the constants part of detModel
  void XercesBuilder::buildConstants(){
    unsigned int i,j;

    if (m_docClient->getConstants() != DOM_Node()){
      Constants* ConstantsBranch = new Constants();
      
      DOM_NodeList childs = m_docClient->getConstants().getChildNodes();
      // version primary ?derived 
      for(i=0;i<childs.getLength();i++){
	if(childs.item(i).getNodeType() != Comment)
	  {
	    std::string str = 
              std::string(xml::Dom::transToChar(childs.item(i).getNodeName()));
	    if(str == "version"){
	      std::string s1,s2;
	      DOM_NamedNodeMap attr=childs.item(i).getAttributes ();
	      s1=
                std::string
                (xml::Dom::transToChar
                 (attr.getNamedItem(DOMString("major")).getNodeValue()));

	      s2=
                std::string
                (xml::Dom::transToChar
                 (attr.getNamedItem(DOMString("minor")).getNodeValue()));

	      ConstantsBranch->setVersion(s1,s2);
	    }
	    else if(str == "primary" || str == "derived"){
	      DOM_NodeList child = childs.item(i).getChildNodes();
	      if (child.getLength()!=0){
		for(j=0;j<child.getLength();j++){
                  if (child.item(j).getNodeType() != Comment) {
                    ConstCategory* cat= new ConstCategory;
                    std::string s1;
                    DOM_NamedNodeMap attrCat=child.item(j).getAttributes();

                    s1=
                      std::string
                      (xml::Dom::transToChar
                       (attrCat.getNamedItem
                        (DOMString("name")).getNodeValue()));

                    cat->setName(s1);
                    // if the actual category has more than one attribute 
                    // it is derived
                    if(attrCat.getLength()>1){ 

                      std::string s2 =
                        std::string(xml::Dom::transToChar
                                    (attrCat.getNamedItem
                                     (DOMString("save")).getNodeValue()));

                      if (s2=="true") {
                        cat->setSave(true);
                      }
                      cat->setPrimary(false);
                    }
                    DOM_Node over =
                      child.item(j).getFirstChild().getFirstChild();
                    s1 =
                      std::string(xml::Dom::transToChar( over.getNodeValue()));
                    cat->setOverview(s1);
                    child.item(j).removeChild(child.item(j).getFirstChild());
                    //elt is the list of prim/const
                    DOM_NodeList elt=child.item(j).getChildNodes();
                    unsigned int k;
                    for(k=0;k<elt.getLength();k++){
                      if (elt.item(k).getNodeType()!=Comment){
                        cat->addConstant(buildConst(elt.item(k)));
                      }//end if 
                    }//end for iterating over k
                    //Insert a new category in Constast object
                    ConstantsBranch->addConstantCategory(cat);
                  } // end if != Comment
		}//end for prim and derived cateory 
	      }
	    }  
	  }
      }//end for                                           
      currentGdd->setConstants(ConstantsBranch);
      MatCollection* materials = new MatCollection();
      currentGdd->setMaterials(materials);
      currentGdd->buildConstantsMap();
    }//end if
  }


  void XercesBuilder::buildMaterials(){

    if (m_docClient->getMaterials() != DOM_Node()){
      MatCollection* materials = currentGdd->getMaterials();

      DOM_NamedNodeMap attrCol=(m_docClient->getMaterials()).getAttributes();

      materials->setVersion(std::string
                            (xml::Dom::transToChar
                             (attrCol.getNamedItem
                              (DOMString("version")).getNodeValue())));

      materials->setDate(std::string
                         (xml::Dom::transToChar
                          (attrCol.getNamedItem
                           (DOMString("date")).getNodeValue())));

      materials->setAuthor(std::string
                           (xml::Dom::transToChar
                            (attrCol.getNamedItem
                             (DOMString("author")).getNodeValue())));

      DOM_NodeList childs = m_docClient->getMaterials().getChildNodes();
      for(unsigned int i=0;i<childs.getLength();i++){
	if(childs.item(i).getNodeType() != Comment)
	  {
	    std::string str = 
              std::string(xml::Dom::transToChar(childs.item(i).getNodeName()));
	    if(str == "element")
	      {
		materials->addMaterial(buildElement((childs.item(i))));
	      }
	    if(str == "composite")
	      {
		materials->addMaterial(buildComposite((childs.item(i))));
	      }
	  }
      }
      
    }//end if
  }

  Element* XercesBuilder::buildElement(DOM_Node e)
  {
    Element* element = new Element();
    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));
    
    element->setName(std::string
                     (xml::Dom::transToChar(el.getAttribute("name"))));

    if(el.getAttribute("symbol") != "")
      element->setSymbol
        (std::string(xml::Dom::transToChar(el.getAttribute("symbol"))));


    if(el.getAttribute("RGB") != "")
      {
	std::string temp = 
          std::string(xml::Dom::transToChar(el.getAttribute("RGB")));
	
      }

    element->setZ(atof(xml::Dom::transToChar(el.getAttribute("z"))));
    element->setAweight(atof(xml::Dom::transToChar
                             (el.getAttribute("aweight"))));

    if(el.getAttribute("density") != "") {
      element->setDensity(atof(xml::Dom::transToChar
                               (el.getAttribute("density"))));
    }

    return element;
  }


  /// Method that build a composite from the relative DOM node
  Composite* XercesBuilder::buildComposite(DOM_Node e){
    unsigned int i, n = 0;

    /// The childs of the node
    DOM_NodeList childs = e.getChildNodes();

    /// Counts the number of addmaterial; we need this loop to avoid comments
    for(i=0;i<childs.getLength();i++)
      {
	if (childs.item(i).getNodeType()!=Comment)
	  n++;
      }

    /// Create a new composite material with n components
    Composite* comp = new Composite(n);

    /// Cast to an element the node
    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));
    

    /// Set the name of the composite
    comp->setName(std::string(xml::Dom::transToChar(el.getAttribute("name"))));
    /// Set its desnsity
    comp->setDensity(atof(xml::Dom::transToChar(el.getAttribute("density"))));

    /// Start to cylce on the addmaterial elements
    for(i=0;i<childs.getLength();i++)
      {
	/// Avoid comments
	if (childs.item(i).getNodeType()!=Comment)
	  {
	    const DOM_Element elCh = 
              DOM_Element(static_cast<const DOM_Element &>(childs.item(i)));
	    
	    /// Get the name of the material to add to the compostite
	    std::string name = 
              std::string(xml::Dom::transToChar
                          (elCh.getAttribute("material")));
	    
	    /// Prepare a material	    
	    Material* mat;
	    
	    /// Search the material in the MatCollection
	    MatCollection* mats = currentGdd->getMaterials();
	    std::map<std::string, Material*> matMap;
	    typedef std::map<std::string, Material*> M;
	    M::iterator mIt;
	    matMap = mats->getMaterials();

	    /// If the material is not yet defined issue an error and abort
	    mIt = matMap.find(name);
	    if(mIt == matMap.end())
	      detAbort("Error in materials: composite uses a not defined material");
	    else mat = mIt->second; 
	    
	    /// Get the child of addmaterial, i.e. a fractionmass or an natoms
	    DOM_Node chi = elCh.getFirstChild();
	    DOM_Element elChi = DOM_Element(static_cast<DOM_Element &>(chi));

	    std::string type = 
              std::string(xml::Dom::transToChar(chi.getNodeName()));
	    /// Add component using the overloaded addComponent of composite
	    if(type == "fractionmass") {
	      comp->addComponent(mat,atof
                                 (xml::Dom::transToChar
                                  (elChi.getAttribute("fraction"))));
            }
	    else if(type == "natoms") {
	      comp->addComponent
                (mat,
                 (unsigned int)atoi(xml::Dom::transToChar
                                    (elChi.getAttribute("n"))));
	    }
	  }
      }
    return comp;
  }



  void XercesBuilder::buildSections()
  {
    if (m_docClient->getSections() != DOM_Node()){
      //      currentGdd->addSection(buildSection(const_cast<DOM_Node*>(m_docClient->getSections())));
      currentGdd->addSection
        (buildSection((m_docClient->getSections())));
    }
    
    currentGdd->buildVolumeMap();

    currentGdd->ResolveReferences();

    currentGdd->buildBoundingBoxes();

  }



  Section* XercesBuilder::buildSection(DOM_Node e)
  {
    unsigned int i;

    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));
    Section* s = new Section(
			     xml::Dom::getAttribute(el, "name"),
			     xml::Dom::getAttribute(el, "version"),
			     xml::Dom::getAttribute(el, "date"),
			     xml::Dom::getAttribute(el, "author"),
			     xml::Dom::getAttribute(el, "fineChoice"),
			     xml::Dom::getAttribute(el, "topVolume")
			     );
    DOM_NodeList childs = e.getChildNodes();
    for(i=0;i<childs.getLength();i++)
      {
	if (childs.item(i).getNodeType()!=Comment){
	  std::string str = 
            std::string(xml::Dom::transToChar(childs.item(i).getNodeName()));
	  if(str == "box")
	    s->addVolume(buildBox(childs.item(i)));
	  if(str == "tubs")
	    s->addVolume(buildTube(childs.item(i)));
	  else if(str == "choice")
	    {
	      Choice* t = buildChoice(childs.item(i));
	      // s->addChoice(t); 
	      s->addVolume(t);
	    }
	  else if(str == "composition")
	    s->addVolume(buildComposition(childs.item(i) ));
	  else if(str=="stackX"||str=="stackY"||str=="stackZ")
	    s->addVolume(buildStack(childs.item(i) ));
	}
      }
    return s;
  }


  Choice* XercesBuilder::buildChoice(DOM_Node e)
  {
    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));
    Choice* b = new Choice(xml::Dom::getAttribute(el, "name"),
			   xml::Dom::getAttribute(el, "default"));

    //    Manager* Manager = Manager::getPointer();

    DOM_NodeList childs = e.getChildNodes();
  
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
	  b->addCaseName(Mode, Volume);
	}
    }
    return b;
  }


  /// \doit Add parameters
  Box* XercesBuilder::buildBox(DOM_Node e)
  {

    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));
    Box* b = new Box(xml::Dom::getAttribute(el, "name"));

    // b->setUnitLength(xml::Dom::transToChar(el.getAttribute("unitLength")));
    // b->setUnitAngle(xml::Dom::transToChar(el.getAttribute("unitAngle")));
    if(el.getAttribute("X") != "0")
      b->setX(atof(xml::Dom::transToChar(el.getAttribute("X"))));
    if(el.getAttribute("Y") != "0")
      b->setY(atof(xml::Dom::transToChar(el.getAttribute("Y"))));
    if(el.getAttribute("Z") != "0")
      b->setZ(atof(xml::Dom::transToChar(el.getAttribute("Z"))));
    b->setMaterial(xml::Dom::transToChar(el.getAttribute("material")));
    b->setDetectorType(atoi(xml::Dom::transToChar
                            (el.getAttribute("detectorType"))));

    if (std::string((xml::Dom::transToChar
                     (el.getAttribute("sensitive")))) == "posHit")
      b->setSensitive(1);
    else if (std::string((xml::Dom::transToChar
                          (el.getAttribute("sensitive")))) == "intHit")
      b->setSensitive(2);
    else
      b->setSensitive(0);

    if (e.hasChildNodes()){
      DOM_Node child=e.getFirstChild();
      Seg* s = new Seg;

      DOM_NamedNodeMap attributelist = child.getAttributes();
      for(unsigned int k=0;k<attributelist.getLength();k++){
	std::string NameAttr=
	  std::string(xml::Dom::transToChar
                      (attributelist.item(k).getNodeName()) );

	char* ValueAttr=xml::Dom::transToChar
          (attributelist.item(k).getNodeValue() );

	if(NameAttr=="axis")s->setAxis(ValueAttr);
	else if(NameAttr=="reason")s->setReason(ValueAttr);
	else if(NameAttr=="nSeg" && "1")s->setnSeg(atoi(ValueAttr));
      }//end for
      b->setSeg(s); 
    }
    return b;
  }



  /// \doit Add parameters
  Tube* XercesBuilder::buildTube(DOM_Node e)
  {

    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));
    Tube* t = new Tube(xml::Dom::getAttribute(el, "name"));


    // t->setUnitLength(xml::Dom::transToChar(el.getAttribute("unitLength")));
    // t->setUnitAngle(xml::Dom::transToChar(el.getAttribute("unitAngle")));
    if(el.getAttribute("RIn") != "0")
      t->setRin(atof(xml::Dom::transToChar(el.getAttribute("RIn"))));
    if(el.getAttribute("ROut") != "0")
      t->setRout(atof(xml::Dom::transToChar(el.getAttribute("ROut"))));
    if(el.getAttribute("Z") != "0")
      t->setZ(atof(xml::Dom::transToChar(el.getAttribute("Z"))));
    t->setMaterial(xml::Dom::transToChar(el.getAttribute("material")));

    if (std::string((xml::Dom::transToChar
                     (el.getAttribute("sensitive")))) == "true")
      t->setSensitive(1);
    else 
      t->setSensitive(0);
    
    return t;
  }



  Composition* XercesBuilder::buildComposition(DOM_Node e){
    unsigned int i;

    Composition* b = new Composition();
    DOM_NamedNodeMap attributelist=e.getAttributes();

    for (i=0;i<attributelist.getLength();i++){
      if (attributelist.item(i).getNodeType()!=Comment){
	DOM_Node currentAttribute=attributelist.item(i);
	std::string value = 
          std::string(xml::Dom::transToChar(currentAttribute.getNodeValue() ));
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

    DOM_NodeList childs = e.getChildNodes();
  
    for(i=0;i<childs.getLength();i++){
      if (childs.item(i).getNodeType()!=Comment)
	b->addPosition(buildPosition( childs.item(i) ));
      
    };

    return b;
  }


  Stack* XercesBuilder::buildStack(DOM_Node e){

    Stack::axisDir st;
    unsigned int i;

    std::string StackDir = 
      std::string(xml::Dom::transToChar( e.getNodeName() ));
    if(StackDir=="stackX")st=Stack::xDir;
    else if(StackDir=="stackY")st=Stack::yDir;
    else st=Stack::zDir;
    
    Stack* s = new Stack(st);
    
    DOM_NamedNodeMap attributelist=e.getAttributes();
    for (i=0;i<attributelist.getLength();i++){
      if (attributelist.item(i).getNodeType()!=Comment){
	DOM_Node currentAttribute=attributelist.item(i);
   
	std::string value = 
          std::string(xml::Dom::transToChar(currentAttribute.getNodeValue() ));
	std::string attributeName = 
	  std::string(xml::Dom::transToChar( currentAttribute.getNodeName() ));
	
	if( attributeName=="name") s->setName(value);
	else if (attributeName=="origin"){
	  if (value=="atStart")s->setOrigin(Stack::atStart);
	  else s->setOrigin(Stack::atCentre);
	}
	else if (attributeName== "parameters")s->setParameters(value);
      }
    }

    DOM_NodeList childs = e.getChildNodes();
  
  
    for(i=0;i<childs.getLength();i++){
      if ( (childs.item(i)).getNodeType()!=Comment){
	StackedPos* p = buildRelativePosition( childs.item(i));
	p->setAxisDir((Stack::axisDir)st);
	s->addPosition(p);
      }
    };
  
    return s;
  }



  SinglePos* XercesBuilder::buildPosition(DOM_Node e){

    std::string posType = std::string(xml::Dom::transToChar(e.getNodeName()));
  
    if (posType=="posXYZ"){
    
      PosXYZ* pos=new PosXYZ();
      setAttributePosition(pos, e);
       
      DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));

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




  StackedPos* XercesBuilder::buildRelativePosition(DOM_Node e){

    std::string posType = std::string(xml::Dom::transToChar(e.getNodeName()));
    //AxisPos,AxisMpos
    if (posType=="axisPos"){
      //   GddaxisPos* pos=new GddaxisPos();
      //   setAttributeRelativePosition(pos, e);
       
      AxisMPos* pos=new AxisMPos();

      setAttributeRelativePosition(pos, e);
       
      DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));

      if(el.getAttribute("shift") != "0")
	pos->setShift0(atof(xml::Dom::transToChar(el.getAttribute("shift"))));
      if(el.getAttribute("gap") != "0")
	pos->setGap0(atof(xml::Dom::transToChar(el.getAttribute("gap"))));

      pos->setNcopy(1);
    

      return pos;
      
    }
  
    else{
      AxisMPos* pos=new AxisMPos();

      setAttributeRelativePosition(pos, e);
       
      DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));

      if(el.getAttribute("shift0") != "0")
	pos->setShift0(atof(xml::Dom::transToChar(el.getAttribute("shift0"))));
      if(el.getAttribute("gap0") != "0")
	pos->setGap0(atof(xml::Dom::transToChar(el.getAttribute("gap0"))));
      if(el.getAttribute("ncopy") != "0")
	pos->setNcopy(atof(xml::Dom::transToChar(el.getAttribute("ncopy"))));
    
      return pos;
    
    }
  }


  void XercesBuilder::setAttributePosition(SinglePos* pos, DOM_Node e)
  {
    DOM_NamedNodeMap attributelistPos=e.getAttributes();

    DOM_Node node;

    node = attributelistPos.getNamedItem(DOMString("volume"));
    pos->setVolumeRef(xml::Dom::transToChar(node.getNodeValue()));

    node = attributelistPos.getNamedItem(DOMString("xrot"));
    pos->setXRot(atof(xml::Dom::transToChar(node.getNodeValue())));

    node = attributelistPos.getNamedItem(DOMString("yrot"));
    pos->setYRot(atof(xml::Dom::transToChar(node.getNodeValue())));

    node = attributelistPos.getNamedItem(DOMString("zrot"));
    pos->setZRot(atof(xml::Dom::transToChar(node.getNodeValue())));

    //    node = attributelistPos.getNamedItem(DOMString("S"));
    //    pos->setS(atof(xml::Dom::transToChar(node.getNodeValue())));

    //  node = attributelistPos.getNamedItem(DOMString("unitLength"));
    // pos->setUnitLength(xml::Dom::transToChar(node.getNodeValue()));

    // node = attributelistPos.getNamedItem(DOMString("unitAngle"));
    // pos->setUnitAngle(xml::Dom::transToChar(node.getNodeValue()));


    if (e.hasChildNodes()){
      DOM_NodeList child=e.getChildNodes();
      for(unsigned int k=0; k<child.getLength() ; k++)
	{
	  if(child.item(k).getNodeType() != Comment) 
	    {
	      IdField* field = new IdField;

	      DOM_NamedNodeMap attributelist = child.item(k).getAttributes();
	      for(unsigned int i=0;i<attributelist.getLength();i++){
	      
		std::string NameAttr=
		  std::string(xml::Dom::transToChar
                              (attributelist.item(i).getNodeName()) );

		char* ValueAttr=xml::Dom::transToChar
                  (attributelist.item(i).getNodeValue() );
		if(NameAttr=="name")field->setName(ValueAttr);
		else if(NameAttr=="step")field->setStep(atof(ValueAttr));
		else if(NameAttr=="value")field->setValue(atof(ValueAttr));
	      }//end for
	      pos->addIdField(field);
	    }//end if
	}
    }
  }

  void XercesBuilder::setAttributeRelativePosition(StackedPos* pos, DOM_Node e)
  {
    DOM_NamedNodeMap attributelistPos=e.getAttributes();

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

    if (e.hasChildNodes()){
      DOM_NodeList child=e.getChildNodes();
      for(unsigned int k=0; k<child.getLength() ; k++)
	{
	  if(child.item(k).getNodeType() != Comment) 
	    {
	      IdField* field = new IdField;

	      DOM_NamedNodeMap attributelist = child.item(k).getAttributes();
	      for(unsigned int i=0;i<attributelist.getLength();i++){
	      
		std::string NameAttr=
		  std::string(xml::Dom::transToChar
                              (attributelist.item(i).getNodeName()) );

		char* ValueAttr=xml::Dom::transToChar
                  ( attributelist.item(i).getNodeValue() );

		if(NameAttr=="name")field->setName(ValueAttr);
		else if(NameAttr=="step")field->setStep(atof(ValueAttr));
		else if(NameAttr=="value")field->setValue(atof(ValueAttr));
	      }//end for
	      pos->addIdField(field);
	    }//end if
	}
    }
 
  }

}


