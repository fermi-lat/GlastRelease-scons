#include <string>

#include "detModel/GDDpurge.h"
#include "detModel/GDD.h"
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
  typedef map<string, GDDvolume*> M;
  M::const_iterator i;
  typeVolume voltype;
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
    voltype=actualVolume->getVolumeType();

    switch(voltype){
    case shape:
      break;
    case composition:
      {
	///number of positions of the current composition
	int numPosition;
	int k;
	   
	/// Static cast of volume to composition
	GDDcomposition* actualComposition=
	  static_cast<GDDcomposition*>(actualVolume);
	GDDvolume* refEnvelope=
	  getVolumeByName(actualComposition->getEnvelopeRef());
	   
	/// Check that the envelope is a simple shape
	if (refEnvelope->getVolumeType() ==shape)
	  {
	    actualComposition->setEnvelope(refEnvelope);
	    numPosition=actualComposition->getPositions().size();
	    GDDanyPosition* currentPos;
	    for(k=0;k<numPosition;k++)
	      {
		currentPos=actualComposition->getPositions()[k];
		GDDvolume* refvol=getVolumeByName(currentPos->getVolumeRef());
		if (refvol!=0){
		  currentPos->setVolume(refvol);
		}
	      }
	  }
	else{
	  //error: envelope must be a shape! We have to put something here 
	}
	break;
      };
    case stack:
      {
	///number of positions of the current stack
	int numPosition;
	int k;
	   
	/// Static cast of volume to stack
	GDDstack* actualStack = static_cast<GDDstack*>(actualVolume);
	   
	numPosition=actualStack->getPositions().size();
	GDDanyRelativePosition* currentRelPos;
	for(k=0;k<numPosition;k++){
	  currentRelPos=actualStack->getPositions()[k];
	  GDDvolume* refvol=getVolumeByName(currentRelPos->getVolumeRef());
	  if (refvol!=0){
	    currentRelPos->setVolume(refvol);
	  }
	}//end for
	break;
      };
    case logical:
      {
	//logicalResolution
	break;
      };
    }//end swich
  }//end for
}

/// This method gives back the total number of volumes in all the sections 
int GDD::getVolumesNumber()
{
  int n=0;
  typedef map<string, GDDvolume*> M;
  M::const_iterator i; 

  for(i=volumeMap.begin();i!=volumeMap.end();i++)
    n++;
  return n;
}

/// This method gives back the total number of sections  
int GDD::getSectionsNumber()
{
  return sections.size();
}

/* This method gives back a GDDchoice* given a name
 * If it does not exist, it returns a null pointer
 */
GDDchoice* GDD::getChoiceByName(string cname)
{
  typedef map<string, GDDchoice*> M;
  M::const_iterator i; 

  i = choiceMap.find(cname);
  if(i == choiceMap.end()) return 0;
  else return i->second;
}

/* This method gives back a GDDvolume* given a name
 * If it does not exist, it returns a null pointer
 */
GDDvolume* GDD::getVolumeByName(string vname)
{
  typedef map<string, GDDvolume*> M;
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
	string mode = GDDmanager::getPointer()->getMode();
	return getVolumeByName(choice->getVolumeNameByMode(mode));
      }
    else return 0;
  }
  else return i->second;
}

/* This method build the volume map and it should be automatically 
 * called by the builder 
 */
void GDD::buildVolumeMap()
{
  unsigned int j;
  int i;
  typedef map<string, GDDvolume*> M;

  for(i=0;i<getSectionsNumber();i++)
    for(j=0;j<sections[i]->getVolumes().size();j++)
      volumeMap.insert(M::value_type(sections[i]->getVolumes()[j]->getName(),
				     sections[i]->getVolumes()[j]));

  buildChoiceMap();
}

/* This method build the volume map and it should be automatically 
 * called by the builder 
 */
void GDD::buildChoiceMap()
{
  unsigned int i,j;
  typedef map<string, GDDchoice*> M;

  for(i=0;i<sections.size();i++)
    for(j=0;j<sections[i]->getChoices().size();j++)
      choiceMap.insert(M::value_type(sections[i]->getChoices()[j]->getName(),
				     sections[i]->getChoices()[j]));
}

/* This method gives back a double given a name
 * If it does not exist, it returns a null pointer
 */
double GDD::getConstantByName(string cname)
{
  typedef map<string, double> M;
  M::const_iterator i; 

  i = constMap.find(cname);
  if(i == constMap.end()) return 0;
  else return i->second;
}

/* This method gives back a string given a name
 * If it does not exist, it returns a null pointer
 */
string GDD::getMaterialByConstantName(string cname)
{
  typedef map<string, string> M;
  M::const_iterator i; 
  
  i = materialMap.find(cname);
  if(i == materialMap.end()) return 0;
  else return i->second;
}

/* This method build the constant map and  the material map,it should be automatically 
 * called by the builder 
 */
void GDD::buildConstantsMap()
{
  unsigned int k,j;

  typedef map<string, double> M1;
  typedef map<string, string> M2;

  for(k=0;k<constants->getConstantCategories().size();k++){
    GDDconstCategory* tmp=constants->getConstantCategories()[k];
    for(j=0;j<tmp->getConsts().size();j++){
      GDDconst* temp = tmp->getConsts()[j];

      switch(temp->getConstType()){
      case i:{
	GDDintConst* tmp1=static_cast<GDDintConst*>(temp);
	double val=(double)tmp1->getValue();
	constMap.insert(M1::value_type(temp->getName(),val));
	break;
      }
      case f:{
	GDDfloatConst* tmp1=static_cast<GDDfloatConst*>(temp);
	double val=(double)tmp1->getValue();
	constMap.insert(M1::value_type(temp->getName(),val));
	break;
      }
      case d:{
	GDDdoubleConst* tmp1=static_cast<GDDdoubleConst*>(temp);
	double val=(double)tmp1->getValue();
	constMap.insert(M1::value_type(temp->getName(),val));
	break;
      }
      case s:{
	GDDstringConst* tmp1=static_cast<GDDstringConst*>(temp);
	materialMap.insert(M2::value_type(temp->getName(),tmp1->getValue()));
	break;
      }
      }//end switch
    }
  }
}

/* This is the recursive Accept function for the sections visitor
 * This method start the visit on the GDD and than recursively
 * start the visit of all the sections it contains
 */
void GDD::Accept(GDDvisitor* v){
  unsigned int i;

  switch(v->getType()){
  case sectionsVisitor:
    {
      GDDsectionsVisitor* vsec = static_cast<GDDsectionsVisitor*>(v);
      vsec->visitGDD(this);
      for(i=0; i<sections.size();i++){
	sections[i]->Accept(vsec);
      }
      break;
    }
  case constantsVisitor:
    //    v->visitGDD(this);
    // constants->Accept(v);
    break;
  }
}

/* This is the non recursive Accept function for the sections visitor
 * This method only start the visit on the GDD and leave to the visitor
 * the option to recurse on its contained elements
 */
void GDD::AcceptNotRec(GDDvisitor* v){
  v->visitGDD(this);
}


void GDD::visitMap(){
  typedef map<string,string>M;
  M::const_iterator i;
  std::cout<<"Material List"<<std::endl;
  for (i=materialMap.begin(); i!=materialMap.end(); i++){
    std::cout<<"constant "<<i->first<<" value: "<<i->second<<std::endl;
  };
  typedef map<string,double>M1;
  M1::const_iterator j;
  std::cout<<"Numeric constants list"<<std::endl;
  for (j=constMap.begin(); j!=constMap.end(); j++){
    std::cout<<"constant "<<j->first<<" value: "<<j->second<<std::endl;
  }
}







