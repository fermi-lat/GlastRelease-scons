// $Header$
#include <vector>

#include "MediaVisitor.h"

#include "detModel/Utilities/Global.h"
#include "detModel/Gdd.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Materials/Element.h"
#include "detModel/Materials/Composite.h"
#include <algorithm> // for accumulate

MediaVisitor::MediaVisitor(IMedia& media):m_media(media)
{
  setRecursive(1);
}

void MediaVisitor::visitElement(detModel::Element* element)
{
  std::string name = element->getName();
  std::string symbol = element->getSymbol();
  double z = element->getZ();
  double a =  element->getAweight();
  double density = element->getDensity(); 

  IMedia::MediaType type;
  std::vector<double> params;
  params.push_back(z);
  params.push_back(a); 

  if (symbol != "") 
    {
      type = IMedia::Element;
    }
  else
    {
      type = IMedia::Material;
      params.push_back(density); 
    }

  m_media.addMaterial(name, type, params, symbol);
  m_materials.push_back(name);
}


void MediaVisitor::visitComposite(detModel::Composite* composite)
{
  std::string name = composite->getName();
  double density = composite->getDensity();
  int ncomponents = composite->getNComponents();

  std::vector <std::string> components;
  std::vector <double> qty(ncomponents);
  IMedia::CompositeType type;

  // Decide the type (fraction or natoms) and prepare the qty vector
  if(composite->isFractions())
    {
      type = IMedia::Fract;
      std::vector <double> fractions = composite->getFractions();
      std::copy(fractions.begin(),fractions.end(),qty.begin());
    }
  else
    {
      type = IMedia::Natoms;
      std::vector <int> natoms = composite->getAtoms();
      std::copy(natoms.begin(),natoms.end(),qty.begin());
    }

  // See if all the components have been visited, otherwise visit them and
  // prepare the components names vector
  std::vector <detModel::Material*> matComponents = composite->getComponents();
  std::vector <detModel::Material*>::iterator m;
  for (m = matComponents.begin(); m != matComponents.end(); m++)
    {
      std::vector<std::string>::iterator where;
      where = std::find(m_materials.begin(),m_materials.end(),(*m)->getName());
      if (where == m_materials.end())
	{
	  (*m)->Accept(this);
	}
      components.push_back((*m)->getName());
    }
  
  m_media.addComposite(name, type, density, components, qty);
  m_materials.push_back(name);
}
