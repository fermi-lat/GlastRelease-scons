#include <string>
#include <iostream>
#include <fstream>
#include <vector>

int main() {

  std::ifstream f("end_seed.dat");

  int nEvents = 10;
  std::string seeds[nEvents];
 
  for(int i = 0; i != nEvents; ++i) {

    string temp;
    getline(f, temp);

    getline(f, seeds[i]);
  }

  if(seeds[0] != seeds[7] || seeds[1] != seeds[5] || seeds[2] != seeds[9] ||
     seeds[3] != seeds[8] || seeds[4] != seeds[6]) {
    std::cout << "test RandomSvc failed!" << std::endl;
    return 1;
  }
  else {
    std::cout << "test RandomSvc ok!" << std::endl;
    return 0;
  }

}


