
#ifndef __PropertiesCheckAlg_H
#define __PropertiesCheckAlg_H 1


#include <GaudiKernel/Algorithm.h>
#include <GaudiKernel/Property.h>
#include <map>

class PropertiesCheckAlg : public Algorithm {
    
public:

    PropertiesCheckAlg( const std::string &, ISvcLocator * ) ;
    
    StatusCode initialize() ;
    StatusCode execute() { return StatusCode::SUCCESS ; }
    StatusCode finalize() ;
    
private: 

    //! clients whose properties must be checked
    StringArrayProperty m_exclude ;
    
    //! count the used properties
    void readPropertyCallBack( Property & ) ;
    
    //! the map of checked properties
    struct PropertyInfo
     {
      std::string client ;
      int used ;
     } ;
    std::map<std::string,PropertyInfo> m_properties ;
} ;

#endif


