#ifndef CHOICE_H
#define CHOICE_H
#include <string>
#include <map>
#include "detModel/Sections/Volume.h"

namespace detModel{

  class Choice :public Volume{
  public:
    Choice(std::string pName, std::string pDef):Volume(pName),defaultName(pDef){;};

    std::string getVolumeNameByMode(std::string);
    void addCaseName(std::string, std::string);

    Volume* getVolumeByMode(std::string);
    void addCase(std::string, Volume*);

    Volume* getDefaultVol(){return defaultVol;}
    void setDefaultVol(Volume* p){defaultVol = p;}

    std::string getDefaultName(){return defaultName;}
    void setDefaultName(std::string p){defaultName = p;}

    std::map <std::string, std::string>* getCasesName(){return &casesName;};

    std::string getMode(){return actualMode;}
    void setMode(std::string mode){actualMode = mode;}

    virtual void Accept(SectionsVisitor*);
    virtual void AcceptNotRec(SectionsVisitor*);

    virtual void buildBB();
    
  private:
    std::string defaultName;
    Volume* defaultVol;

    std::string actualMode;

    std::map <std::string, std::string> casesName;
    std::map <std::string, Volume*> cases;
  };
}
#endif //CHOICE_H
