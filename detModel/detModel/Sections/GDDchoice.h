#ifndef GDDCHOICE_H
#define GDDCHOICE_H
#include <map>
#include "detModel/Sections/GDDvolume.h"

class GDDchoice {
    public:
        GDDchoice(string pname):name(pname){};
        string getName(){return name;};
        string getVolumeNameByMode(string mode);
        void addCase(string mode, string volume);

    private:
        string name;
        map <string, string> cases;
};
#endif //GDDCHOICE_H
