#include "detModel/Materials/Composite.h"
#include "detModel/Materials/Material.h"
#include "detModel/Management/MaterialsVisitor.h"

namespace detModel{

  void Composite::Accept(MaterialsVisitor* v){
    
    v->visitComposite(this);
    
  }



  void Composite::AcceptNotRec(MaterialsVisitor* v){
    v->visitComposite(this);
  }
  
  
  unsigned int Composite::getAtom(std::string mName){
    if (isAtoms())
      {
	unsigned int i;
	
	for(i=0;i<nComponents;i++)
	  if (components[i]->getName() == mName)
	    return atoms[i];
      }
    return 0;
  }


  double Composite::getFraction(std::string mName){
    if (isFractions())
      {
	unsigned int i;
	
	for(i=0;i<nComponents;i++)
	  if (components[i]->getName() == mName)
	    return fractions[i];
      }
    return 0.0;    
  }
  
  void Composite::addComponent(Material* m, unsigned int n){

    if (fractions.size())
      {    
	detAbort("Attention: you are mixing fractions and atoms number in a composite material");
      }
    
    if(components.size() < nComponents)
      {
	components.push_back(m);
	atoms.push_back(n);
	
      }
    else
      {
	detAbort("Attention: you are adding too many components to a composite");
      }

  }

  void Composite::addComponent(Material* m, double d){

    if (atoms.size())
      {    
	detAbort("Attention: you are mixing fractions and atoms number in a composite material");
      }
    
    if(components.size() < nComponents)
      {
	if((d<0.0) || (d>1.0))
	  {
	    detAbort("Attention:  fractions in composite materials must be comprised between 0.0 and 1.0");
	  }

	components.push_back(m);
	fractions.push_back(d);
	
      }
    else
      {
	detAbort("Attention: you are adding too many components to a composite");
      }

    if(components.size() == nComponents)
      {
	double temp = 0.0;
	unsigned int i;

	for(i=0;i<nComponents;i++)
	  temp += fractions[i];

	if (temp != 1.0)
	  detAbort("Attention:  the fractions of a composite material does not add to 1.0");
      }

  }

}
