#ifndef GDDCHOICE_H
#define GDDCHOICE_H
#include <map>
#include "detModel/Sections/GDDvolume.h"

class GDDchoice {
    public:
        GDDchoice(std::string pname):name(pname){};
        std::string getName(){return name;};
        std::string getVolumeNameByMode(std::string mode);
        void addCase(std::string mode, std::string volume);

    private:
        std::string name;
        std::map <std::string, std::string> cases;
};
#endif //GDDCHOICE_H
