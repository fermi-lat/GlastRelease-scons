#include <string>
#include <iostream>

#include "detModel/Gdd.h"
#include "detModel/Sections/Ensemble.h"
#include "detModel/Utilities/Purge.h"
#include "detModel/Sections/Shape.h"
#include "detModel/Sections/Section.h"
#include "detModel/Sections/Choice.h"
#include "detModel/Sections/Volume.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Management/Manager.h"
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Management/ConstantsVisitor.h"
#include "detModel/Constants/Constants.h"
#include "detModel/Constants/IntConst.h"
#include "detModel/Constants/FloatConst.h"
#include "detModel/Constants/DoubleConst.h"
#include "detModel/Constants/StringConst.h"
#include "detModel/Constants/ConstCategory.h"


namespace detModel{

  /// This is the destructor of Gdd 
  Gdd::~Gdd(){
    if (constants)delete constants;
    /// With this call we recursively delete the vector of sections
    purge(sections);
  }

  /// This method resolve the volumes references in the XML file
  void Gdd::ResolveReferences(){
    typedef std::map<std::string, Volume*> M;
    M::const_iterator i;
    Volume* actualVolume;
  
    /// This part resolves the topVolumes
    for(unsigned int j=0;j<sections.size();j++)
      {
	Volume* vol = getVolumeByName(sections[j]->getTopVolumeRef());
	if(vol)
	  sections[j]->setTopVolume(vol);
      }

    /// This part resolves the envelopes and positions volumes
    for(i=volumeMap.begin();i!=volumeMap.end();i++){
      actualVolume=i->second;

      if (Ensemble* actualEns = dynamic_cast<Ensemble*>(actualVolume))
	{
	  ///number of positions of the current composition
	  int numPosition;
	  int k;
	
	  Position* currentPos;
	
	  numPosition=actualEns->getPositions().size();
	  for(k=0;k<numPosition;k++)
	    {
	      currentPos=actualEns->getPositions()[k];
	      Volume* refvol=getVolumeByName(currentPos->getVolumeRef());
	      if (refvol!=0){
		currentPos->setVolume(refvol);
	      }
	    }
	
	  if(Composition* actualComp = dynamic_cast<Composition*>(actualEns))
	    {
	      Shape* refEnvelope=
		getShapeByName(actualComp->getEnvelopeRef());	  
	      if (refEnvelope)
		actualComp->setEnvelope(refEnvelope);	  
	    }
	}
    }
  }

  /// This method gives back the total number of volumes in all the sections 
  int Gdd::getVolumesNumber()
  {
    int n=0;
    typedef std::map<std::string, Volume*> M;
    M::const_iterator i; 

    for(i=volumeMap.begin();i!=volumeMap.end();i++)
      n++;
    return n;
  }

  /// This method gives back the total number of constants  
  int Gdd::getConstantsNumber()
  {
    int n=0;
    typedef std::map<std::string, Const*> M;
    M::const_iterator i; 

    for(i=constMap.begin();i!=constMap.end();i++)
      n++;
    return n;
  }

  Const* Gdd::getOrderedConst(int ii)
  {
    int n=0;
    typedef std::map<std::string, Const*> M;
    M::const_iterator i; 

    for(i=constMap.begin();i!=constMap.end();i++){
      if (n == ii) return i->second;
      else n++;
    }
    return 0;
  }

  Volume* Gdd::getOrderedVolume(int ii)
  {
    int n=0;
    typedef std::map<std::string, Volume*> M;
    M::const_iterator i; 

    for(i=volumeMap.begin();i!=volumeMap.end();i++){
      if (n == ii) return i->second;
      else n++;
    }

    return 0;
  }

  /// This method gives back the total number of sections  
  int Gdd::getSectionsNumber()
  {
    return sections.size();
  }

  /* This method gives back a Choice* given a name
   * If it does not exist, it returns a null pointer
   */
  Choice* Gdd::getChoiceByName(std::string cname)
  {
    typedef std::map<std::string, Choice*> M;
    M::const_iterator i; 

    i = choiceMap.find(cname);
    if(i == choiceMap.end()) return 0;
    else return i->second;
  }

  /* This method gives back a Volume* given a name
   * If it does not exist, it returns a null pointer
   */
  Volume* Gdd::getVolumeByName(std::string vname)
  {
    typedef std::map<std::string, Volume*> M;
    M::const_iterator i; 

    i = volumeMap.find(vname);
    if(i == volumeMap.end()){
      Choice* choice = getChoiceByName(vname);
    
      /* We check that the name correspond to a choice
       * In that case we return a pointer to a Volume
       * depending to the mode, contained in the manager.
       */
      if(choice)
	{
	  std::string mode = Manager::getPointer()->getMode();
	  return getVolumeByName(choice->getVolumeNameByMode(mode));
	}
      else return 0;
    }
    else return i->second;
  }

  /* This method gives back a Volume* given a name
   * If it does not exist, it returns a null pointer
   */
  Shape* Gdd::getShapeByName(std::string vname)
  {
  
    if( getVolumeByName(vname))
      return dynamic_cast<Shape*>(getVolumeByName(vname));
    else return 0;
  }

  /* This method build the volume map and it should be automatically 
   * called by the builder 
   */
  void Gdd::buildVolumeMap()
  {
    unsigned int j;
    int i;
    typedef std::map<std::string, Volume*> M;

    for(i=0;i<getSectionsNumber();i++)
      for(j=0;j<sections[i]->getVolumes().size();j++)
	volumeMap.insert(M::value_type(sections[i]->getVolumes()[j]->getName(),
				       sections[i]->getVolumes()[j]));
    buildChoiceMap();
  }

  void Gdd::buildBoundingBoxes()
  {
    unsigned int j;
    int i;

    for(i=0;i<getSectionsNumber();i++)
      for(j=0;j<sections[i]->getVolumes().size();j++)
	sections[i]->getVolumes()[j]->buildBB();
  }

  /* This method build the choices map and it is called by the buildVolumeMap
   */
  void Gdd::buildChoiceMap()
  {
    unsigned int i,j;
    typedef std::map<std::string, Choice*> M;

    for(i=0;i<sections.size();i++)
      for(j=0;j<sections[i]->getChoices().size();j++)
	choiceMap.insert(M::value_type(sections[i]->getChoices()[j]->getName(),
				       sections[i]->getChoices()[j]));

    /// We fill the modeNames vector
    for(j=0;j<sections[0]->getChoices().size();j++)
      modeNames.push_back(sections[0]->getChoices()[j]->getName());

  }

  ConstCategory* Gdd::getConstCategoryByName(std::string cname){

    std::vector<ConstCategory*>::const_iterator i;
  
    Constants* con = getConstants();
    std::vector <ConstCategory*> vc = con->getCategories();
    for(i=vc.begin();i!=vc.end();i++)
      if ((*i)->getName() == cname) return (*i); 

    return 0;
  }


  /* This method gives back a Gddconst* given a name
   * If it does not exist, it returns a null pointer
   */
  Const* Gdd::getConstByName(std::string cname)
  {
    typedef std::map<std::string, Const*> M;
    M::const_iterator i; 

    i = constMap.find(cname);
    if(i == constMap.end()) return 0;
    else return i->second;
  }

  /* This method gives back a double given a name
   * If it does not exist, it returns a null pointer
   */
  double Gdd::getNumConstByName(std::string cname)
  {
    typedef std::map<std::string, double> M;
    M::const_iterator i; 

    i = constNumMap.find(cname);
    if(i == constNumMap.end()) return 0;
    else return i->second;
  }

  /* This method gives back a string given a name
   * If it does not exist, it returns a null pointer
   */
  std::string Gdd::getCharConstByName(std::string cname)
  {
    typedef std::map<std::string, std::string> M;
    M::const_iterator i; 
  
    i = constCharMap.find(cname);
    if(i == constCharMap.end()) return 0;
    else return i->second;
  }

  /* This method build the constant map and  the material map,it should be automatically 
   * called by the builder 
   */
  void Gdd::buildConstantsMap()
  {
    unsigned int k,j,m;

    typedef std::map<std::string, double> M1;
    typedef std::map<std::string, std::string> M2;
    typedef std::map<std::string, Const *> M3;

    for(k=0;k<constants->getCategories().size();k++){
      ConstCategory* tmp=constants->getCategories()[k];
      for(j=0;j<tmp->getConsts().size();j++){
	Const* temp = tmp->getConsts()[j];

	constMap.insert(M3::value_type(temp->getName(),temp));

	if (IntConst* tmp1 = dynamic_cast<IntConst*>(temp)){ 
	  double val=(double)tmp1->getValue();
	  constNumMap.insert(M1::value_type(temp->getName(),val));
	}
	else if (FloatConst* tmp1 = dynamic_cast<FloatConst*>(temp)){ 
	  double val=(double)tmp1->getValue();
	  constNumMap.insert(M1::value_type(temp->getName(),val));
	}
	else if (DoubleConst* tmp1 = dynamic_cast<DoubleConst*>(temp)){       
	  double val=(double)tmp1->getValue();
	  constNumMap.insert(M1::value_type(temp->getName(),val));
	}
	else if (StringConst* tmp1 = dynamic_cast<StringConst*>(temp)){       
	  constCharMap.insert(M2::value_type(temp->getName(),tmp1->getValue()));
	  /// Here we fill the materialNames vector
	  if (temp->getConstMeaning() == Const::mat)
	    {
	      bool flag=1;
	      for(m=0;m<materialNames.size();m++)
		if (materialNames[m]==tmp1->getValue()) flag = 0;
	      if (flag) materialNames.push_back(tmp1->getValue());
	    }
	}
      }
    }
  }

  /* This is the recursive Accept function for the sections visitor
   * This method start the visit on the Gdd and than recursively
   * start the visit of all the sections it contains
   */
  void Gdd::Accept(Visitor* v){
    unsigned int i;

    if(SectionsVisitor* vsec = dynamic_cast<SectionsVisitor*>(v))
      {
	vsec->visitGdd(this);
	for(i=0; i<sections.size();i++){
	  sections[i]->Accept(vsec);
	}
      }
    else if(ConstantsVisitor* vcon = dynamic_cast<ConstantsVisitor*>(v))
      {
	vcon->visitGdd(this);
	constants->Accept(vcon);
      }
  }


  /* This is the non recursive Accept function for the sections visitor
   * This method only start the visit on the Gdd and leave to the visitor
   * the option to recurse on its contained elements
   */
  void Gdd::AcceptNotRec(Visitor* v){
    v->visitGdd(this);
  }

}








