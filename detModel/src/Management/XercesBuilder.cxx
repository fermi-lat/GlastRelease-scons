// #include "xmlUtil/Substitute.h"
// #include "xmlUtil/Arith.h"
#include "xmlUtil/docMan/GDDDocMan.h"
#include "xmlUtil/id/IdDict.h"

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>
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

  XERCES_CPP_NAMESPACE_USE

  XercesBuilder::XercesBuilder()
  {
    m_docClient = new DMDocClient();
    
    xmlUtil::GDDDocMan* pGDDMan = xmlUtil::GDDDocMan::getPointer();
    pGDDMan->regClient("constants", m_docClient);
    pGDDMan->regClient("section", m_docClient);    
    pGDDMan->regClient("materials", m_docClient);    
    pGDDMan->regClient("idDict", m_docClient);    
  }
  
  XercesBuilder::~XercesBuilder()
  {
    delete m_docClient;
  }

  bool XercesBuilder::parseFile(std::string nameFile){
//    unsigned int iSec;
    xmlUtil::GDDDocMan* pGDDMan = xmlUtil::GDDDocMan::getPointer();
    
    bool parseOk = pGDDMan->parse(nameFile);
    
    /// We start detModel stuff retrieving the manager
    Manager* man = Manager::getPointer();

    /// We point to the current Gdd
    currentGdd = man->getGdd();
  
    /// Set some info on the parsed file
    if (parseOk) {
      currentGdd->setCVSid(pGDDMan->getCVSid());
      currentGdd->setDTDversion(pGDDMan->getDTDversion());
    }
    return parseOk;
  }  


  /*
   *  This method build a constant
   */
  Const* XercesBuilder::buildConst(DOMElement* e){
    using xml::Dom;

    std::string name;
    std::string typeOfConst;
    std::string ut;

    //    DOM_NamedNodeMap attr = e.getAttributes();

    name = Dom::getAttribute(e, "name");
    ut = Dom::getAttribute(e, "uType");
    std::string elementName = Dom::getTagName(e);

    if (elementName=="prim"){
      typeOfConst=Dom::getAttribute(e,"type");

      Const* c;
      try {
        if (typeOfConst=="int"){
          c=new IntConst;
          int val= Dom::getIntAttribute(e, "value");
          ((IntConst*)c)->setValue(val);
        }
        else if (typeOfConst=="float"){
          c=new FloatConst;
          double val = Dom::getDoubleAttribute(e, "value");
          float  valFloat = val;
          ((FloatConst*)c)->setValue(valFloat);
        }
        else if (typeOfConst=="double"){
          c=new DoubleConst;
          double val=
            Dom::getDoubleAttribute(e, "value");
          ((DoubleConst*)c)->setValue(val);
        }
        else if (typeOfConst=="string"){
          c=new StringConst;
          std::string val=Dom::getAttribute(e, "value");
          ((StringConst*)c)->setValue(val);
        }
        else return 0;
      }  // end try
      catch (xml::DomException ex) {
        std::cerr << "from detModel::XercesBuild::BuildConst" << std::endl
                  << ex.getMsg() << std::endl;
        throw (ex);
      }

      c->setName(name);
      c->setConstMeaning(ut);
      try {
        c->setNote(Dom::getText(e->getFirstChild()));
      }
      catch (xml::DomException ex) {
        // this isn't terribly fatal, so don't rethrow
        std::cerr << "From detModel::XercesBuilder::getConst " << std::endl
                  << ex.getMsg() << std::endl;
      }
      return c;
    }//end if prim
    else {
      //e is a const element.  Can be any of the numeric types
      // For backward compatibility, if no "type" attribute, treat
      // it as a double
      std::string valType;
      Const* c;

      try {
        if (!Dom::hasAttribute(e, "type")) valType = std::string("double");
        else valType = Dom::getAttribute(e, "type");
        if (valType == "int") {
          c = new IntConst;
          int val = Dom::getIntAttribute(e, "value");
          ((IntConst*)c)->setValue(val);
        }
        else if (valType == "float") {
          c = new FloatConst;
          double valDouble = Dom::getDoubleAttribute(e, "value");
          float val = valDouble;
          ((FloatConst*)c)->setValue(val);
        }
        else if (valType == "double") {  
          c = new DoubleConst;
          double val = Dom::getDoubleAttribute(e, "value");
          ((DoubleConst*)c)->setValue(val);
        }  
        else return 0;  // or, better, maybe throw exception?
      }                     // end try
      catch (xml::DomException ex) {
        std::cerr << "From detModel::XercesBuilder::buildConst" << std::endl
                  << ex.getMsg() << std::endl;
        throw ex;
      }

      c->setName(name);
      c->setConstMeaning(ut);
    
      DOMElement* notes = Dom::findFirstChildByName(e, "notes");
      if (notes != 0) {
        c->setNote(Dom::getText(notes->getFirstChild()));
      }
      return c;
    }//end else (const, not prim)
  }

  /// This methods build the constants part of detModel
  void XercesBuilder::buildConstants(){
    using xml::Dom;

    if (m_docClient->getConstants() != 0){
      Constants* constantsBranch = new Constants();

      const DOMNode* node = m_docClient->getConstants();
      const DOMElement* constants = 
        static_cast<const DOMElement*>(node);
      // Content model for <constants> is (version, primary, derived?)

      {
        DOMElement* version = Dom::findFirstChildByName(constants, "version");
        std::string majorVersion = Dom::getAttribute(version, "major");
        std::string minorVersion = Dom::getAttribute(version, "minor");
        constantsBranch->setVersion(majorVersion, minorVersion);
      }
      
      buildCategories(constantsBranch, 
                      Dom::findFirstChildByName(constants, "primary"), true);
      buildCategories(constantsBranch, 
                      Dom::findFirstChildByName(constants, "derived"), false);

      currentGdd->setConstants(constantsBranch);
      MatCollection* materials = new MatCollection();
      currentGdd->setMaterials(materials);
      currentGdd->buildConstantsMap();
    }//end if
  }

  void XercesBuilder::buildCategories(Constants* constantsBranch, 
                                      DOMElement* parent, bool primary) {
    std::vector<DOMElement*> cats;

    using xml::Dom;

    // Could be called with null arg if no derived
    if (parent == 0) return;  
    
    Dom::getChildrenByTagName(parent, "*", cats);
    unsigned nCat = cats.size();

    for (unsigned iCat = 0; iCat < nCat; iCat++) {  // loop over categories
      ConstCategory* cat = new ConstCategory;
      cat->setName(Dom::getAttribute(cats[iCat], "name"));
      cat->setPrimary(primary);

      if (!primary) {
        std::string save = Dom::getAttribute(cats[iCat], "save");
        cat->setSave(save == "true");
      }

      DOMElement* overview = Dom::findFirstChildByName(cats[iCat], "overview");
      cat->setOverview(Dom::getText(overview->getFirstChild()));

      // loop over remaining children, the constants
      DOMElement* constElt = Dom::getSiblingElement(overview);
      while (constElt != 0) {
        cat->addConstant(buildConst(constElt));
        constElt = Dom::getSiblingElement(constElt);
      }

      constantsBranch->addConstantCategory(cat);
    }
  } 

  void XercesBuilder::buildMaterials(){
    using xml::Dom;

    if (m_docClient->getMaterials() == 0) return;

    MatCollection* materials = currentGdd->getMaterials();

    const DOMNode* node = m_docClient->getMaterials();
    const DOMElement* matsElt = 
      static_cast<const DOMElement*>(node);
      
    materials->setVersion(Dom::getAttribute(matsElt, "version"));
    materials->setDate(Dom::getAttribute(matsElt, "date"));
    materials->setAuthor(Dom::getAttribute(matsElt, "author"));


    std::vector<DOMElement*> children;
    Dom::getChildrenByTagName(matsElt, "*", children);

    unsigned nChild = children.size();
    for (unsigned iChild = 0; iChild < nChild; iChild++) {
      std::string tagName = Dom::getTagName(children[iChild]);
      if (tagName == "element") {
        materials->addMaterial(buildElement(children[iChild]) );
      }
      else if (tagName == "composite") {
        materials->addMaterial(buildComposite(children[iChild]) );
      }
    }
  }

  void XercesBuilder::buildIdDictionary() {
    if (m_docClient->getIdDictionary() != 0){
      DOMNode*    dictNode = m_docClient->getIdDictionary();
      DOMElement* idDictElt = 
	static_cast<DOMElement*>(dictNode);

      xmlUtil::IdDict *dict = new xmlUtil::IdDict(idDictElt);

      if (dict->isValid()) {
	currentGdd->setIdDictionary(dict);
      }
    }
  }

  Element* XercesBuilder::buildElement(DOMElement* e)
  {
    using xml::Dom;

    Element* element = new Element();
    
    element->setName(Dom::getAttribute(e, "name"));

    std::string symbol = Dom::getAttribute(e, "symbol");
    if (symbol.size() > 0) element->setSymbol(symbol);
    try {
      element->setZ(Dom::getDoubleAttribute(e, ("z")));
      element->setAweight(Dom::getDoubleAttribute(e, ("aweight")));
      if (Dom::hasAttribute(e, "density")) {
        element->setDensity(Dom::getDoubleAttribute(e, ("density")));
      }
    }
    catch (xml::DomException ex) {
      std::cerr << "From detModel::XercesBuilder::buildElement" << std::endl
                << ex.getMsg() << std::endl;
      throw ex;
    }
    return element;
  }


  /// Method that build a composite from the relative DOM node
  Composite* XercesBuilder::buildComposite(DOMElement* e){
    using xml::Dom;

    std::vector<DOMElement*> addMats;

    Dom::getChildrenByTagName(e, "addmaterial", addMats);
    unsigned int n = addMats.size();

    /// Create a new composite material with n components
    Composite* comp = new Composite(n);

    /// Set the name of the composite
    comp->setName(Dom::getAttribute(e, "name"));

    /// Set its desnsity
    try {
      comp->setDensity(Dom::getDoubleAttribute(e, "density"));
    }
    catch (xml::DomException ex) {
      std::cerr << "From detModel::XercesBuilder::buildComposite" 
                << std::endl << ex.getMsg() << std::endl;
      throw ex;
    }

    // Start to cylce on the addmaterial elements
    for (unsigned int i = 0; i < n; i++) {
      std::string name = Dom::getAttribute(addMats[i], "material");
      Material* mat;

      // Search the material in the MatCollection
      MatCollection* mats = currentGdd->getMaterials();
      std::map<std::string, Material*> matMap;
      typedef std::map<std::string, Material*> M;
      M::iterator mIt;
      matMap = mats->getMaterials();
      
      // If the material is not yet defined issue an error and abort
      mIt = matMap.find(name);
      if(mIt == matMap.end())
        detAbort("Error in materials: composite uses an undefined material");
      else mat = mIt->second; 

      DOMElement* child = Dom::getFirstChildElement(addMats[i]);
      std::string type = Dom::getTagName(child);
      
      try {
        if (type == "fractionmass") {
          comp->addComponent(mat, 
                             Dom::getDoubleAttribute(child, "fraction"));
        }
        else if (type == "natoms") {
          comp->addComponent(mat, 
                             (unsigned int) Dom::getIntAttribute(child, "n"));
        }
      }
      catch (xml::DomException ex) {
        std::cerr << "From detModel::XercesBuilder::buildComposite" 
                  << std::endl << ex.getMsg() << std::endl;
        throw ex;
      }

    }
    return comp;
  }
  


  void XercesBuilder::buildSections()   {
    if (m_docClient->getSections() != 0){
      currentGdd->addSection
        (buildSection((m_docClient->getSections())));
    }
    
    currentGdd->buildVolumeMap();
    
    currentGdd->ResolveReferences();
    
    currentGdd->buildBoundingBoxes();
    
  }
  
  Section* XercesBuilder::buildSection(const DOMNode* e)
  {
    using xml::Dom;

    const DOMElement* el = static_cast<const DOMElement *>(e);
    Section* s = new Section(
			     Dom::getAttribute(el, "name"),
			     Dom::getAttribute(el, "version"),
			     Dom::getAttribute(el, "date"),
			     Dom::getAttribute(el, "author"),
			     Dom::getAttribute(el, "fineChoice"),
			     Dom::getAttribute(el, "topVolume")
			     );
    std::vector<DOMElement*> children;

    Dom::getChildrenByTagName(el, "*", children);
    
    unsigned int n = children.size();

    for (unsigned int i=0; i < n; i++) {
      std::string tagName = Dom::getTagName(children[i]);

      if (tagName == "box") s->addVolume(buildBox(children[i]));
      else if (tagName == "tubs") s->addVolume(buildTube(children[i]));
      else if (tagName == "choice") s->addVolume(buildChoice(children[i]));
      else if (tagName == "composition") {
        s->addVolume(buildComposition(children[i]));
      }
      else if ((tagName == "stackX")  ||
               (tagName == "stackY")  ||
               (tagName == "stackZ")   ) {
        s->addVolume(buildStack(children[i]));
      }
    }
    return s;
  }


  Choice* XercesBuilder::buildChoice(DOMElement* e)
  {
    using xml::Dom;

    //    DOM_Element el = DOM_Element(static_cast<DOM_Element &>(e));
    Choice* b = new Choice(Dom::getAttribute(e, "name"),
			   Dom::getAttribute(e, "default"));

    //    Manager* Manager = Manager::getPointer();

    //    DOM_NodeList childs = e.getChildNodes();
    std::vector<DOMElement*> children;
    Dom::getChildrenByTagName(e, "*", children);
    unsigned int n = children.size();

    for (unsigned int i=0; i < n; i++) {
	  
      std::string Mode = Dom::getAttribute(children[i], "mode");

      std::string Volume = Dom::getAttribute(children[i], "volume");

      b->addCaseName(Mode, Volume);
    }
    return b;
  }


  /// \doit Add parameters
  Box* XercesBuilder::buildBox(DOMElement* e)
  {
    using xml::Dom;
    Box* b = new Box(Dom::getAttribute(e, "name"));

    std::string Xatt = Dom::getAttribute(e, "X"); 
    std::string Yatt = Dom::getAttribute(e, "Y"); 
    std::string Zatt = Dom::getAttribute(e, "Z"); 
    try {
      if (Dom::hasAttribute(e, "X") ) {
        b->setX(Dom::getDoubleAttribute(e, "X"));
      }
      if (Dom::hasAttribute(e, "Y") ) {
        b->setY(Dom::getDoubleAttribute(e, "Y"));
      }
      if (Dom::hasAttribute(e, "Z") ) {
        b->setZ(Dom::getDoubleAttribute(e, "Z"));
      }
    }
    catch (xml::DomException ex) {
      std::cerr << "From detModel::XercesBuilder::buildBox  " << std::endl
                  << ex.getMsg() << std::endl;
        throw ex;
    }


    b->setMaterial(Dom::getAttribute(e, "material"));
    b->setDetectorType(Dom::getIntAttribute(e, "detectorType"));

    std::string sensAtt = Dom::getAttribute(e, "sensitive");
    if (sensAtt == "posHit") b->setSensitive(1);
    else if (sensAtt == "intHit") b->setSensitive(2);
    else b->setSensitive(0);

    // A box may be segmented
    DOMElement* segElt = Dom::getFirstChildElement(e);
    if (segElt != 0) {
      Seg* s = new Seg;

      // axis and reason attributes are required by dtd
      s->setAxis(Dom::getAttribute(segElt, "axis"));
      s->setReason(Dom::getAttribute(segElt, "reason"));

      // nSeg is defaulted; nSegREF shouldn't be present
      // by the time this code gets to the xml source
      try {
        s->setnSeg(Dom::getIntAttribute(segElt, "nSeg"));
      }
      catch (xml::DomException ex) {
        std::cerr << "From detModel::XercesBuilder::buildBox  " << std::endl
                  << ex.getMsg() << std::endl;
        throw ex;
      }
                 
      b->setSeg(s); 
    }
    return b;
  }



  /// \doit Add parameters
  Tube* XercesBuilder::buildTube(DOMElement* e)
  {
    using xml::Dom;

    Tube* t = new Tube(Dom::getAttribute(e, "name"));


    try {
      double attVal = Dom::getDoubleAttribute(e, "RIn");
      if (attVal != 0) t->setRin(attVal);
      attVal = Dom::getDoubleAttribute(e, "ROut");
      if (attVal != 0) t->setRout(attVal);
      attVal = Dom::getDoubleAttribute(e, "Z");
      if (attVal != 0) t->setZ(attVal);
    }
    catch (xml::DomException ex) {
      std::cerr << "From detModel::XercesBuilder::buildTube" << std::endl
                << ex.getMsg() << std::endl;
      throw ex;
    }

    t->setMaterial(Dom::getAttribute(e, "material"));

    std::string sensAtt = Dom::getAttribute(e, "sensitive");
    if (sensAtt == "posHit") t->setSensitive(1);
    else if (sensAtt == "intHit") t->setSensitive(2);
    else t->setSensitive(0);
    return t;
  }


    
  Composition* XercesBuilder::buildComposition(DOMElement* e)  {
    using xml::Dom;

    Composition* b = new Composition();

    // Deal with the attributes
    b->setName(Dom::getAttribute(e, "name"));
    b->setEnvelopeRef(Dom::getAttribute(e, "envelope"));
    if (Dom::hasAttribute(e, "parameters")) {
      b->setParameters(Dom::getAttribute(e, "parameters"));
    }
     
    std::vector<DOMElement*> posElts;
    Dom::getChildrenByTagName(e, "*", posElts);
    unsigned int n = posElts.size();
    
    for (unsigned int i = 0; i < n; i++) {
      b->addPosition(buildPosition(posElts[i]));
    }
  
    return b;
  }


  Stack* XercesBuilder::buildStack(DOMElement* e) {

    using xml::Dom;

    Stack::axisDir st;

    std::string StackDir = Dom::getTagName(e);

    if (StackDir=="stackX") st=Stack::xDir;
    else if (StackDir=="stackY") st=Stack::yDir;
    else st=Stack::zDir;
    
    Stack* s = new Stack(st);

    // name and origin attributes are always there
    s->setName(Dom::getAttribute(e, "name"));
    std::string origin = Dom::getAttribute(e, "origin");
    if (origin == "atStart") s->setOrigin(Stack::atStart);
    else s->setOrigin(Stack::atCentre);

    if (Dom::hasAttribute(e, "parameters")) {
      s->setParameters(Dom::getAttribute(e, "parameters"));
    }

    std::vector<DOMElement*> posElts;
    
    Dom::getChildrenByTagName(e, "*", posElts);

    unsigned int n = posElts.size();

    for (unsigned int i = 0; i < n; i++) {
      StackedPos* p = buildRelativePosition(posElts[i]);
      p->setAxisDir((Stack::axisDir)st);
      s->addPosition(p);
    }
    return s;
  }



  SinglePos* XercesBuilder::buildPosition(DOMElement* e) {
    using xml::Dom;

    std::string posType =  Dom::getTagName(e);
  
    if (posType=="posXYZ"){
    
      PosXYZ* pos=new PosXYZ();
      setAttributePosition(pos, e); 

      try {
        if (Dom::hasAttribute(e, "X")) {
          pos->setX(Dom::getDoubleAttribute(e, "X"));
        }
        if (Dom::hasAttribute(e, "Y")) {
          pos->setY(Dom::getDoubleAttribute(e, "Y"));
        }
        if (Dom::hasAttribute(e, "Z")) {
          pos->setZ(Dom::getDoubleAttribute(e, "Z"));
        }
      }
      catch (xml::DomException ex) {
        std::cerr << "From detModel::XercesBuilder::buildPosition" << std::endl
                  << ex.getMsg() << std::endl;
        throw ex;
      }

      return pos;
    }
    else if (posType=="posRPhiZ"){
      /// \todo Add posRPhiZ
      return 0;
    }
    else  return 0; // default
  }

  StackedPos* XercesBuilder::buildRelativePosition(DOMElement* e) {

    using xml::Dom;
    std::string posType = Dom::getTagName(e);

    //AxisPos,AxisMpos
    if (posType=="axisPos"){
       
      AxisMPos* pos=new AxisMPos();

      setAttributeRelativePosition(pos, e);

      if (Dom::hasAttribute(e, "shift")) {
        pos->setShift0(Dom::getDoubleAttribute(e, "shift"));
      }
      if (Dom::hasAttribute(e, "gap")) {
        pos->setGap0(Dom::getDoubleAttribute(e, "gap"));
      }

      pos->setNcopy(1);
    
      return pos;
    }
  
    else{
      AxisMPos* pos=new AxisMPos();

      setAttributeRelativePosition(pos, e);
       
      try {
        if (Dom::hasAttribute(e, "shift0")) {
          pos->setShift0(Dom::getDoubleAttribute(e, "shift0"));
        }
        if (Dom::hasAttribute(e, "gap0")) {
          pos->setGap0(Dom::getDoubleAttribute(e, "gap0"));
        }
        if (Dom::hasAttribute(e, "ncopy")) {
          pos->setNcopy(Dom::getDoubleAttribute(e, "ncopy"));
        }
      }
      catch (xml::DomException ex) {
        std::cerr << "From detModel::XercesBuilder::buildRelativePosition"
                  << std::endl << ex.getMsg() << std::endl;
      }
      return pos;
    }
  }


  void XercesBuilder::setAttributePosition(SinglePos* pos, DOMElement* e)
  {
    using xml::Dom;

    pos->setVolumeRef(Dom::getAttribute(e, "volume"));

    try {
      pos->setXRot(Dom::getDoubleAttribute(e, "xrot"));
      pos->setYRot(Dom::getDoubleAttribute(e, "yrot"));
      pos->setZRot(Dom::getDoubleAttribute(e, "zrot"));
    }
    catch (xml::DomException ex) {
      std::cerr << "From detModel::XercesBuilder::setAttributePosition"
                << std::endl << ex.getMsg() << std::endl;
    }
                      

    // Handle id field children, if any
    setIdFields(pos, e);

  }

  void XercesBuilder::setAttributeRelativePosition(StackedPos* pos, 
                                                   DOMElement* e)  {
    using xml::Dom;
    pos->setVolumeRef(Dom::getAttribute(e, "volume"));

    try {
      pos->setDx(Dom::getDoubleAttribute(e, "dX"));
      pos->setDy(Dom::getDoubleAttribute(e, "dY"));
      pos->setDz(Dom::getDoubleAttribute(e, "dZ"));
      pos->setRotation(Dom::getDoubleAttribute(e, "rotation"));
      pos->setGap(Dom::getDoubleAttribute(e, "gap"));
      pos->setShift(Dom::getDoubleAttribute(e, "shift"));
    }
    catch (xml::DomException ex) {
      std::cerr << "From detModel::XercesBuild::setAttributeRelativePosition"
                << std::endl << ex.getMsg() << std::endl;
      throw ex;
    }
    // Handle id field children, if any
    setIdFields(pos, e);

  }     // end setAttributeRelativePosition

  void XercesBuilder::setIdFields(Position* pos, DOMElement* e) {
    using xml::Dom;

    std::vector<DOMElement*> fields;
    Dom::getChildrenByTagName(e, "idField", fields);

    unsigned int n = fields.size();

    for (unsigned int i = 0; i < n; i++) {
      IdField* field = new IdField;
      
      DOMElement* curField = fields[i];
      field->setName(Dom::getAttribute(curField, "name"));
      try {
        field->setStep(Dom::getIntAttribute(curField, "step"));
        field->setValue(Dom::getIntAttribute(curField, "value"));
      }
      catch (xml::DomException ex) {
        std::cerr << "From detModel::XercesBuild::setIdFields: " 
                  << std::endl << ex.getMsg() << std::endl;
        throw ex;
      }
      pos->addIdField(field);
    }

  }
} // end namespace detModel



