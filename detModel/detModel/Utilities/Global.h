#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>

namespace detModel{

  inline void detAbort(std::string s){
    std::cout << s << std::endl;
    std::cout << "#################" << std::endl;
    std::cout << "Aborting detModel" << std::endl;
    exit(0);
  }

}
#endif // GLOBAL_H

