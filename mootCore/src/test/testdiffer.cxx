#include <iostream>
#include <cstdio>
#include <string>

bool differ(const std::string& f1, const std::string& f2);

int main(void) {

  std::string* f1 = new
    std::string("/scratch/glastsoft/MOOT/v1r0/MOOD/xml/mood_test.xml");
  std::string* f2 = new
    std::string("/scratch/glastsoft/MOOT/v1r0/MOOD/xml/mood_test.xml~");
  
  if (differ(*f1, *f2)) {
    std::cout << "files f1,f2 differed" << std::endl;
  }
  else {
    std::cout << "files f1,f2 were identical" << std::endl;
  }

  if (differ(*f1, *f1)) {
    std::cout << "files f1,f1 differed" << std::endl;
  }
  else {
    std::cout << "files f1,f1 were identical" << std::endl;
  }

}
bool differ(const std::string& f1, const std::string& f2) {
  std::string cmd("cmp -s ");
  cmd += f1;
  cmd += std::string(" ");
  cmd += f2;
  int ret = system(cmd.c_str());
  return (ret != 0);
}  


