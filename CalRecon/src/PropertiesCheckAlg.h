
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
    StringArrayProperty m_targets ;
    
    //! count the used properties
    void readPropertyCallBack( Property & ) ;
    
    //! the map of used properties
    std::map<std::string,int> m_properties ;
} ;

#endif


