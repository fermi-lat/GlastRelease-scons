#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Gdd.h"
#include "detModel/Sections/Section.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Tube.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/Ensemble.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Sections/IdField.h"
#include "detModel/Sections/Seg.h"
#include "detModel/Management/CountMaterial.h"
#include "detModel/Management/Manager.h"
#include "detModel/Utilities/Global.h"

namespace detModel{

  CountMaterial::CountMaterial(std::string nvol)
  {
    setRecursive(0);
    m_actualVolume = nvol;
  };

  void CountMaterial::visitGdd(Gdd* gdd)
  {
    typedef std::vector<Section*>sec;
    std::vector <Section*>::iterator i;

    /// We iterate on all the sections (actually there should be only 1)
    sec s = gdd->getSections();
    for(i=s.begin(); i!=s.end();i++)
      (*i)->AcceptNotRec(this);

    printTable();
  }

  void  CountMaterial::visitSection(Section* section)
  {
    Volume* vol=0;
    std::vector <Volume*>::iterator v;

    /// If no volume has been specified in the constructor we search for the TopVolume
    if (m_actualVolume == "")
      {
	std::vector<Volume*> volumes = section->getVolumes();
      
	for(v=volumes.begin(); v!=volumes.end(); v++){
	  if(Ensemble* ens = dynamic_cast<Ensemble*>(*v))
	    {
	      if (ens == section->getTopVolume())
		vol = ens;
	    }
	}
      }
    else /// otherwise we use the volume speicified in the constructor
      {
	Manager* manager = Manager::getPointer();
	if (manager->getGdd()->getVolumeByName(m_actualVolume))
	  vol = manager->getGdd()->getVolumeByName(m_actualVolume);
	else /// we abort if such a volume does not exists
	  detAbort("No such a volume!!");
      }

    /// Here we start the visit of the hierarchy.
    vol->AcceptNotRec(this);
  }

  void  CountMaterial::visitEnsemble(Ensemble* ensemble)
  {
    std::vector <Position*>::iterator i;
    typedef std::vector<Position*> pos;
  
    /// Here the positioned volumes are visited
    pos p = ensemble->getPositions();
    for(i=p.begin(); i!=p.end();i++)
      (*i)->AcceptNotRec(this);
      
  }
  
  void CountMaterial::visitPosXYZ(PosXYZ* pos)
  {
    pos->getVolume()->Accept(this);
  }

  void CountMaterial::visitAxisMPos(AxisMPos* pos)
  {
    unsigned int i;
    for(i=0;i<pos->getNcopy();i++)
      pos->getVolume()->Accept(this);
  }
  
  void CountMaterial::visitBox(Box* b)
  {
    m_matNumber[b->getMaterial()]++;
  }
  
  void CountMaterial::visitTube(Tube* t)
  {
    m_matNumber[t->getMaterial()]++;
  }

  void CountMaterial::printTable()
  {
    std::cout << "Materials usage summary: " << std::endl;
    for( std::map<std::string, double>::const_iterator it = m_matNumber.begin(); 
	 it!= m_matNumber.end(); ++it)
      {
	std::cout << "Material " 
		  << (*it).first << "\t" 
		  << (*it).second  
		  << std::endl;
      }
  }

}













