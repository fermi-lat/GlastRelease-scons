#include <vector>

#include "G4MaterialsVisitor.h"

#include "G4Material.hh"
#include "globals.hh"

#include "detModel/Utilities/Global.h"
#include "detModel/Gdd.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Materials/Element.h"
#include "detModel/Materials/Composite.h"


G4MaterialsVisitor::G4MaterialsVisitor()
{
  setRecursive(1);
}

void G4MaterialsVisitor::visitGdd(detModel::Gdd*)
{
}

void G4MaterialsVisitor::visitMatCollection(detModel::MatCollection*)
{
}

void G4MaterialsVisitor::visitElement(detModel::Element* element)
{
  G4String name, symbol;
  G4double a, z, density;
  G4Element* ptElement = 0;
  G4Material* ptMaterial = 0;
  
  if (((ptElement = G4Element::GetElement((G4String) element->getName())) ||
	(ptMaterial = G4Material::GetMaterial((G4String) element->getName()))))  
    return;
  
  name = (G4String) element->getName();
  symbol = (G4String) element->getSymbol();
  z = (G4double) element->getZ();
  a = (G4double) element->getAweight()*g/mole;
  density = (G4double) element->getDensity()*g/cm3;

  if (symbol != "") {
    G4Element* el = new G4Element(name, symbol, z, a);
  }
  else {
      if(name==std::string("Vacuum"))return; // right thing?
      if( z==0  ) { z=6; a=12*g/mole;density=2.265*g/cm3;} // kluge: make it into carbon

    G4Material* el = new G4Material(name, z, a, density);
  }
}


void G4MaterialsVisitor::visitComposite(detModel::Composite* composite)
{
  G4String name;
  G4double density;
  G4int ncomponents;
  G4Material* ptMaterial = 0;
  
  if ((ptMaterial = G4Material::GetMaterial((G4String) composite->getName())))  
    return;  

  name = (G4String) composite->getName();
  density = (G4double) composite->getDensity()*g/cm3;
  ncomponents = (G4int) composite->getNComponents();

  G4Material* mat = new G4Material(name, density, ncomponents);

  std::vector <detModel::Material*> matComponents = composite->getComponents();
  std::vector <detModel::Material*>::iterator m;

  std::vector <int> natoms = composite->getAtoms();
  std::vector <double> fractions = composite->getFractions();
  int i = 0;    //index of the natoms and fractions vectors
 
  for (m = matComponents.begin(); m != matComponents.end(); m++)
    {
      if ( detModel::Element* el = dynamic_cast<detModel::Element*>(*m))
	{
	  G4Element* ptElement = 0;
	  G4Material* ptMaterial = 0;
	  
	  if (!((ptElement = G4Element::GetElement((G4String) el->getName())) ||
		(ptMaterial = G4Material::GetMaterial((G4String) el->getName()))))
	    {
	      el->Accept(this);
	      if (!(ptElement = G4Element::GetElement((G4String) el->getName())))
		{
		  ptMaterial = G4Material::GetMaterial((G4String) el->getName());
		}
	    }
	  if (ptElement)
	    {
	      if (composite->isFractions())
		{
		  mat->AddElement(ptElement, (G4double) fractions[i]);
		}
	      else
		{
		  mat->AddElement(ptElement, (G4int) natoms[i]);
		}
	      i++;
	    }
	  else if (ptMaterial)
	    {
	      if (composite->isFractions())
		{
		  mat->AddMaterial(ptMaterial, (G4double) fractions[i]);
		}
	      else
		{
		  mat->AddMaterial(ptMaterial, (G4int) natoms[i]);
		}
	      i++;
	    }
	  else
	    {
	      detModel::detAbort("Error: material not defined");
	    }
	}
      else
	{
	  detModel::Composite* comp = dynamic_cast<detModel::Composite*>(*m);
	  G4Material* ptMaterial = G4Material::GetMaterial(comp->getName());
	  if (!ptMaterial)
	    {
	      comp->Accept(this);
	      ptMaterial = G4Material::GetMaterial((G4String) comp->getName());
	    }
	  if (ptMaterial)
	    {
	      mat->AddMaterial(ptMaterial, (G4double) fractions[i]);
	      i++;
	    }
	  else
	    {
	      detModel::detAbort("Error: material not defined");
	    }
	}
    }
}
		   
	 
	   
	      
  

  
  


