// $Header $ 

#ifndef IDDICTMAN_H
#define IDDICTMAN_H

#include <vector>
#include <string>
namespace xmlUtil {
  class IdDict;

  class IdDictMan {
  public:
    enum RetCode {
      FAIL = 0,
      SUCCESS = 1,
      NOTFOUND = 2,
      DUPLICATE = 3
    };
    static IdDictMan * IdDictMan::getPointer();

    IdDict * findDict(std::string name);

    RetCode  registerDict(IdDict *dict);

    RetCode  removeDict(std::string name);

  protected:
    IdDictMan() {}

    std::vector<IdDict *> dictDir;
    static IdDictMan *me;

  };
}
#endif
