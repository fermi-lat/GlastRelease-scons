#include "G4Media.h"

#include "G4Material.hh"
#include "globals.hh"

#include <numeric> // for accumulate


void G4Media::addMaterial(std::string name, 
			  MediaType type,
			  const DoubleVector& params, std::string symbol)
{
  if (type == IMedia::Element)
    {
      new G4Element(name, symbol, params[0], params[1]*g/mole);    
   }
  else 
    {// This is temporary waiting for a new material category with pressure and temp
      if(name==std::string("Vacuum"))
	{
	  new G4Material(name, params[0], params[1]*g/mole, params[2]*g/cm3, kStateGas,
			 2.73*kelvin,3.e-18*pascal);
	}
      else
	new G4Material(name, params[0], params[1]*g/mole, params[2]*g/cm3);
 
    }
}


void G4Media::addComposite(std::string name, 
			   CompositeType type,
			   double density,
			   const StringVector& components, 
			   const DoubleVector& qty)
{
  G4Element* ptElement = 0;
  G4Material* ptMaterial = 0;
  
  G4Material* mat = new G4Material(name, density*g/cm3, components.size());
  
  for(unsigned int i=0;i<components.size();i++)
    {
      if (!(ptElement = G4Element::GetElement(components[i])))
	{
	  ptMaterial = G4Material::GetMaterial(components[i]);
	}
      
      if (ptElement)
	{
	  if (type == IMedia::Fract)
	    {
	      double tot = std::accumulate(qty.begin(),qty.end(), 0.0);
	      mat->AddElement(ptElement, (G4double) qty[i]/tot);
	    }
	  else
	    {
	      mat->AddElement(ptElement, (G4int) qty[i]);
	    }
	}
      else if (ptMaterial)
	{
	  if (type == IMedia::Fract)
	    {
	      double tot = std::accumulate(qty.begin(),qty.end(), 0.0);
	      mat->AddMaterial(ptMaterial, (G4double) qty[i]/tot);
	    }
	  else
	    {
	      mat->AddMaterial(ptMaterial, (G4int) qty[i]);
	    }
	}
    }
}
