// $Header $
 
#ifndef CALIBUTIL_CLIENTOBJECT_H
#define CALIBUTIL_CLIENTOBJECT_H

namespace calibUtil {

  typedef struct stowerRC { unsigned int row; unsigned int col;} towerRC;
  enum uniL {TOP, BOT};

  class ClientObject {
  public: 
    
    /// Get data
    virtual unsigned int readData(towerRC towerId, unsigned int trayNum, 
                      uniL uniLayer, std::vector<unsigned int> v) = 0;

  };

}/// end of namespace calibUtil

#endif


