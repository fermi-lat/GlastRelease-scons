#include <map>

#include "detModel/Materials/Material.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Management/MaterialsVisitor.h"
#include "detModel/Utilities/Color.h"
#include "detModel/Utilities/ColorCreator.h"

namespace detModel{


  void MatCollection::addMaterialName(std::string name)
  {
    typedef std::map<std::string,Material*>M;
    typedef std::map<std::string,Color*>C;

    materials.insert(M::value_type(name,0));  
    materialColors.insert(C::value_type(name,new Color()));  
  }


  void MatCollection::addMaterial(Material* m)
  {
    typedef std::map<std::string,Material*>M;
    M::iterator i; 

    i = materials.find(m->getName());
    if(i == materials.end())
      materials.insert(M::value_type(m->getName(),m));
    else i->second = m; 
  }

  void MatCollection::Accept(MaterialsVisitor* v){
    typedef std::map<std::string,Material*>M;
    M::const_iterator i; 
    
    v->visitMatCollection(this);
    
    for(i=materials.begin(); i != materials.end();i++){
      if(i->second)
	{
	  (i->second)->Accept(v);
	}
      else std::cout << "Warning: " << i->first << " has not been defined!" << std::endl;
    }
  }

  void MatCollection::AcceptNotRec(MaterialsVisitor* v){
    v->visitMatCollection(this);
  }


  void MatCollection::generateColor()
  {
    typedef std::map<std::string,Color*>M;
    M::iterator i; 
    unsigned int n = 0;

    ColorCreator* cColor = new ColorCreator(getMaterialsNumber());
    
    for(i=materialColors.begin();i!=materialColors.end();i++)
      {
	(i->second)->setRed(cColor->getColor(n)->getRed());
	(i->second)->setGreen(cColor->getColor(n)->getGreen());
	(i->second)->setBlue(cColor->getColor(n)->getBlue());
	(i->second)->setTra(cColor->getColor(n)->getTra());
	n++;
      }
  }

  Color* MatCollection::getMaterialColor(std::string name)
  {
    typedef std::map<std::string,Color*>M;
    M::iterator i; 

    i = materialColors.find(name);
    if(i == materials.end())
      return 0;
    else return i->second;
  }

  void MatCollection::setMaterialColor(std::string name, float r, float g, float b, float t)
  {
    /// Retrive an iterator to navigate in the material colors map
    typedef std::map<std::string,Color*>M;
    M::iterator i;
    
    i = materialColors.find(name);
    /// Search for name matching and set color if found
    if(i == materials.end())
      return 0;
    else (i->second)->setColor(r,g,b,t);
  }

  void MatCollection::setMaterialColor(std::string name, Color* c)
  {
    /// Retrive an iterator to navigate in the material colors map
    typedef std::map<std::string,Color*>M;
    M::iterator i;
    
    i = materialColors.find(name);
    /// Search for name matching and set color if found
    if(i == materials.end())
      return 0;
    else (i->second)->setColor(c->getRed(),c->getGreen(),c->getBlue(),c->getTra());
  }


  void MatCollection::setMaterialTransparency(std::string name,float t)
  {
    /// Retrive an iterator to navigate in the material colors map
    typedef std::map<std::string,Color*>M;
    M::iterator i;
    
    i = materialColors.find(name);
    /// Search for name matching and set transparency if found
    if(i == materials.end())
      return 0;
    else (i->second)->setTra(t);
  }

  void MatCollection::setAllMaterialTransparency(float t)
  {
    /// Retrive an iterator to navigate in the material colors map
    typedef std::map<std::string,Color*>M;
    M::iterator i;
    
    /// Set transparency for all materials
    for(i=materialColors.begin();i!=materialColors.end();i++)
      (i->second)->setTra(t);      
  }

}


