#include <string>
#include <iostream>

#include "detModel/GDD.h"
#include "detModel/Sections/GDDensamble.h"
#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Sections/GDDshape.h"
#include "detModel/Sections/GDDsection.h"
#include "detModel/Sections/GDDchoice.h"
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Management/GDDconstantsVisitor.h"
#include "detModel/Constants/GDDconstants.h"
#include "detModel/Constants/GDDintConst.h"
#include "detModel/Constants/GDDfloatConst.h"
#include "detModel/Constants/GDDdoubleConst.h"
#include "detModel/Constants/GDDstringConst.h"
#include "detModel/Constants/GDDconstCategory.h"

/// This is the destructor of GDD 
GDD::~GDD(){
  if (constants)delete constants;
  /// With this call we recursively delete the vector of sections
  purge(sections);
}

/// This method resolve the volumes references in the XML file
void GDD::ResolveReferences(){
  typedef std::map<std::string, GDDvolume*> M;
  M::const_iterator i;
  GDDvolume* actualVolume;
  
  /// This part resolves the topVolumes
  for(unsigned int j=0;j<sections.size();j++)
    {
      GDDvolume* vol = getVolumeByName(sections[j]->getTopVolumeRef());
      if(vol)
	sections[j]->setTopVolume(vol);
    }

  /// This part resolves the envelopes and positions volumes
  for(i=volumeMap.begin();i!=volumeMap.end();i++){
    actualVolume=i->second;

    if (GDDensamble* actualEns = dynamic_cast<GDDensamble*>(actualVolume))
      {
	///number of positions of the current composition
	int numPosition;
	int k;
	
	GDDposition* currentPos;
	
	numPosition=actualEns->getPositions().size();
	for(k=0;k<numPosition;k++)
	  {
	    currentPos=actualEns->getPositions()[k];
	    GDDvolume* refvol=getVolumeByName(currentPos->getVolumeRef());
	    if (refvol!=0){
	      currentPos->setVolume(refvol);
	    }
	  }
	
	if(GDDcomposition* actualComp = dynamic_cast<GDDcomposition*>(actualEns))
	  {
	    GDDshape* refEnvelope=
	      getShapeByName(actualComp->getEnvelopeRef());	  
	    if (refEnvelope)
	      actualComp->setEnvelope(refEnvelope);	  
	  }
      }
  }
}

/// This method gives back the total number of volumes in all the sections 
int GDD::getVolumesNumber()
{
  int n=0;
  typedef std::map<std::string, GDDvolume*> M;
  M::const_iterator i; 

  for(i=volumeMap.begin();i!=volumeMap.end();i++)
    n++;
  return n;
}

/// This method gives back the total number of constants  
int GDD::getConstantsNumber()
{
  int n=0;
  typedef std::map<std::string, GDDconst*> M;
  M::const_iterator i; 

  for(i=constMap.begin();i!=constMap.end();i++)
    n++;
  return n;
}

GDDconst* GDD::getOrderedConst(int ii)
{
  int n=0;
  typedef std::map<std::string, GDDconst*> M;
  M::const_iterator i; 

  for(i=constMap.begin();i!=constMap.end();i++){
    if (n == ii) return i->second;
    else n++;
  }
  return 0;
}

GDDvolume* GDD::getOrderedVolume(int ii)
{
  int n=0;
  typedef std::map<std::string, GDDvolume*> M;
  M::const_iterator i; 

  for(i=volumeMap.begin();i!=volumeMap.end();i++){
    if (n == ii) return i->second;
    else n++;
  }

    return 0;
}

/// This method gives back the total number of sections  
int GDD::getSectionsNumber()
{
  return sections.size();
}

/* This method gives back a GDDchoice* given a name
 * If it does not exist, it returns a null pointer
 */
GDDchoice* GDD::getChoiceByName(std::string cname)
{
  typedef std::map<std::string, GDDchoice*> M;
  M::const_iterator i; 

  i = choiceMap.find(cname);
  if(i == choiceMap.end()) return 0;
  else return i->second;
}

/* This method gives back a GDDvolume* given a name
 * If it does not exist, it returns a null pointer
 */
GDDvolume* GDD::getVolumeByName(std::string vname)
{
  typedef std::map<std::string, GDDvolume*> M;
  M::const_iterator i; 

  i = volumeMap.find(vname);
  if(i == volumeMap.end()){
    GDDchoice* choice = getChoiceByName(vname);
    
    /* We check that the name correspond to a choice
     * In that case we return a pointer to a GDDvolume
     * depending to the mode, contained in the manager.
     */
    if(choice)
      {
	std::string mode = GDDmanager::getPointer()->getMode();
	return getVolumeByName(choice->getVolumeNameByMode(mode));
      }
    else return 0;
  }
  else return i->second;
}

/* This method gives back a GDDvolume* given a name
 * If it does not exist, it returns a null pointer
 */
GDDshape* GDD::getShapeByName(std::string vname)
{
  
  if( getVolumeByName(vname))
    return dynamic_cast<GDDshape*>(getVolumeByName(vname));
  else return 0;
}

/* This method build the volume map and it should be automatically 
 * called by the builder 
 */
void GDD::buildVolumeMap()
{
  unsigned int j;
  int i;
  typedef std::map<std::string, GDDvolume*> M;

  for(i=0;i<getSectionsNumber();i++)
    for(j=0;j<sections[i]->getVolumes().size();j++)
      volumeMap.insert(M::value_type(sections[i]->getVolumes()[j]->getName(),
				     sections[i]->getVolumes()[j]));
  buildChoiceMap();
}

void GDD::buildBoundingBoxes()
{
  unsigned int j;
  int i;

  for(i=0;i<getSectionsNumber();i++)
    for(j=0;j<sections[i]->getVolumes().size();j++)
      sections[i]->getVolumes()[j]->buildBB();
}

/* This method build the choices map and it is called by the buildVolumeMap
 */
void GDD::buildChoiceMap()
{
  unsigned int i,j;
  typedef std::map<std::string, GDDchoice*> M;

  for(i=0;i<sections.size();i++)
    for(j=0;j<sections[i]->getChoices().size();j++)
      choiceMap.insert(M::value_type(sections[i]->getChoices()[j]->getName(),
				     sections[i]->getChoices()[j]));

  /// We fill the modeNames vector
  for(j=0;j<sections[0]->getChoices().size();j++)
    modeNames.push_back(sections[0]->getChoices()[j]->getName());

}

GDDconstCategory* GDD::getConstCategoryByName(std::string cname){

  std::vector<GDDconstCategory*>::const_iterator i;
  
  GDDconstants* con = getConstants();
  vector <GDDconstCategory*> vc = con->getCategories();
  for(i=vc.begin();i!=vc.end();i++)
      if ((*i)->getName() == cname) return (*i); 

  return 0;
}


/* This method gives back a GDDconst* given a name
 * If it does not exist, it returns a null pointer
 */
GDDconst* GDD::getConstByName(std::string cname)
{
  typedef std::map<std::string, GDDconst*> M;
  M::const_iterator i; 

  i = constMap.find(cname);
  if(i == constMap.end()) return 0;
  else return i->second;
}

/* This method gives back a double given a name
 * If it does not exist, it returns a null pointer
 */
double GDD::getNumConstByName(std::string cname)
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
std::string GDD::getCharConstByName(std::string cname)
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
void GDD::buildConstantsMap()
{
  unsigned int k,j,m;

  typedef std::map<std::string, double> M1;
  typedef std::map<std::string, std::string> M2;
  typedef std::map<std::string, GDDconst *> M3;

  for(k=0;k<constants->getCategories().size();k++){
    GDDconstCategory* tmp=constants->getCategories()[k];
    for(j=0;j<tmp->getConsts().size();j++){
      GDDconst* temp = tmp->getConsts()[j];

      constMap.insert(M3::value_type(temp->getName(),temp));

      if (GDDintConst* tmp1 = dynamic_cast<GDDintConst*>(temp)){ 
	double val=(double)tmp1->getValue();
	constNumMap.insert(M1::value_type(temp->getName(),val));
      }
      else if (GDDfloatConst* tmp1 = dynamic_cast<GDDfloatConst*>(temp)){ 
	double val=(double)tmp1->getValue();
	constNumMap.insert(M1::value_type(temp->getName(),val));
      }
      else if (GDDdoubleConst* tmp1 = dynamic_cast<GDDdoubleConst*>(temp)){       
	double val=(double)tmp1->getValue();
	constNumMap.insert(M1::value_type(temp->getName(),val));
      }
      else if (GDDstringConst* tmp1 = dynamic_cast<GDDstringConst*>(temp)){       
	constCharMap.insert(M2::value_type(temp->getName(),tmp1->getValue()));
	/// Here we fill the materialNames vector
	if (temp->getConstMeaning() == GDDconst::mat)
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
 * This method start the visit on the GDD and than recursively
 * start the visit of all the sections it contains
 */
void GDD::Accept(GDDvisitor* v){
  unsigned int i;

  if(GDDsectionsVisitor* vsec = dynamic_cast<GDDsectionsVisitor*>(v))
    {
      vsec->visitGDD(this);
      for(i=0; i<sections.size();i++){
	sections[i]->Accept(vsec);
      }
    }
  else if(GDDconstantsVisitor* vcon = dynamic_cast<GDDconstantsVisitor*>(v))
    {
      vcon->visitGDD(this);
      constants->Accept(vcon);
    }
}


/* This is the non recursive Accept function for the sections visitor
 * This method only start the visit on the GDD and leave to the visitor
 * the option to recurse on its contained elements
 */
void GDD::AcceptNotRec(GDDvisitor* v){
  v->visitGDD(this);
}









