#ifndef CHOICE_H
#define CHOICE_H
#include <map>
#include "detModel/Sections/Volume.h"

namespace detModel{

  class Choice {
  public:
    Choice(std::string pname):name(pname){};
    std::string getName(){return name;};
    std::string getVolumeNameByMode(std::string mode);
    void addCase(std::string mode, std::string volume);
    
  private:
    std::string name;
    std::map <std::string, std::string> cases;
  };
}
#endif //CHOICE_H
