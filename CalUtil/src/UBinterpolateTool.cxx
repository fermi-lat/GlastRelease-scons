/**
 *
 * @class UBinterpolateTool
 *
 * @brief interpolation class for the unbiased energy algorithm
 *
 * This class provides the interpolation values for Carmelo's implementation
 * of the unbiased energy algorithm.
 * The function setBiasMap() reads in the interpolation table, the function
 * interpolate() returns an interpolated value for a given logE-zDir pair.
 *
 * @author Michael Kuss, Carmelo Sgro'
 *
 * 
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

#include "CalUtil/IUBinterpolateTool.h"

#include "facilities/Util.h"

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>

class UBinterpolateTool : public AlgTool, virtual public IUBinterpolateTool
{

public:
  /// Standard Gaudi Tool interface constructor
  UBinterpolateTool(const std::string& type,
                const std::string& name,
                const IInterface* parent);
  /// Destructor.
  virtual ~UBinterpolateTool() {}
  /// @brief Initialization of the tool.
  StatusCode initialize();
  /// load map
  void addBiasMap(std::string mapName, std::string calibFileName);
  /// returns an interpolated value for zDir and LogE
  float interpolate( std::string mapName,  float logE,  float zDir) ;

private:
    /// these are map of items, the keys are the id of the energy method name
    /// 
    std::map<std::string,std::string> m_calibFileName;
    /// vector of LogE
    std::map<std::string,std::vector<float> > m_x;
    /// vector of ZDir
    std::map<std::string,std::vector<float> > m_y;
    /// matrix of values
    std::map<std::string,std::vector<float> > m_z;

  float getX(std::string mapName, int i) { 
    if (m_x.find(mapName) != m_x.end() ) { return m_x[mapName][i]; }
    return -2.; }

  float getY( std::string mapName,  int i)  { 
    if (m_y.find(mapName) != m_y.end() ) { return m_y[mapName][i]; }
    return -2.; }
  
  /// returns a value from the interpolation table
  float getZ( std::string mapName,  int i,  int j)  { 
    if (m_z.find(mapName) != m_z.end() ) { return m_z[mapName][i+m_x[mapName].size()*j]; }
    return -2.; }

};

DECLARE_TOOL_FACTORY(UBinterpolateTool);

// Standard Constructor
UBinterpolateTool::UBinterpolateTool(const std::string& type,
				     const std::string& name,
				     const IInterface* parent):
  AlgTool(type, name, parent)
{
  declareInterface<IUBinterpolateTool>(this);
  return;
}

StatusCode UBinterpolateTool::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;
  // Instantiate the message logger.
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "UBinterpolate is initializing..," << endreq;

  return sc;
}

// Set bias map and parse it to fill data vectors
void UBinterpolateTool::addBiasMap(std::string mapName, std::string calibFileName)  {


    m_calibFileName[mapName] = calibFileName;

    facilities::Util::expandEnvVar(&m_calibFileName[mapName]);
    std::ifstream fin(m_calibFileName[mapName].c_str());
    if ( !fin ) {
        std::cerr << "Couldn't open " << m_calibFileName[mapName] << std::endl;
        std::exit(1);
    }
    // tmp vector: fill these and then add to the map
    std::vector<float>  a_x;
    std::vector<float>  a_y;
    std::vector<float>  a_z;
    std::string s;
    int counter = 0;
    while ( !fin.eof() ) {
        std::getline(fin, s);
        if ( s.size() == 0 )
            continue;
        std::stringstream line(s);
        while ( !line.eof() ) {
            float f;
            line >> f;
            switch ( counter ) {
            case 0:
                a_x.push_back(f);
                break;
            case 1:
                a_y.push_back(f);
                break;
            default:
                a_z.push_back(f);
                break;
            }
        }
        ++counter;
    }
    // fill the map 
    m_x[mapName] = a_x;
    m_y[mapName] = a_y;
    m_z[mapName] = a_z;
}


float UBinterpolateTool::interpolate( std::string mapName,  float logE,  float zDir)  {
  
  static  int OUTSIDE = 0;

  // checking the bias map has been added
  if (m_calibFileName.find(mapName) == m_calibFileName.end() ) {
    std::cout << "Map still needs to be added: map: " << mapName << std::endl;
    return OUTSIDE;
  }

    // checking if we are outside of the calibrated area
    if ( logE < getX(mapName, 0) )
        return OUTSIDE;
    if ( logE > getX(mapName, m_x[mapName].size()-1) )
        return OUTSIDE;
    if ( zDir > getY(mapName, 0) )
        return OUTSIDE;
    if ( zDir < getY(mapName, m_y[mapName].size()-1) )
        return OUTSIDE;

    // checking where we are
    unsigned int i = 0;
    unsigned int j = 0;
    for ( i=1; i<m_x[mapName].size()-1; ++i )
        if ( logE < getX(mapName, i) )
            break;
    for ( j=1; j<m_y[mapName].size()-1; ++j )
        if ( zDir > getY(mapName, j) )
            break;

     float x   = logE;
     float y   = zDir;
     float x0  = getX(mapName, i-1);
     float x1  = getX(mapName, i);
     float y0  = getY(mapName, j-1);
     float y1  = getY(mapName, j);
     float z00 = getZ(mapName, i-1,j-1);
     float z01 = getZ(mapName, i-1,j);
     float z10 = getZ(mapName, i  ,j-1);
     float z11 = getZ(mapName, i  ,j);

    // interpolating first in x
     float z0 = ( z00 * ( x1 - x ) - z10 * ( x0 - x ) ) / ( x1 - x0 );
     float z1 = ( z01 * ( x1 - x ) - z11 * ( x0 - x ) ) / ( x1 - x0 );

    // interpolating z0 and z1 in y
     float z = ( z0 * ( y1 - y ) - z1 * ( y0 - y ) ) / ( y1 - y0 );

     /*
    std::cout << "UBinterpolate: map: " << mapName << " file: " 
    << m_calibFileName[mapName] << std::endl;
    std::cout << "we start " << logE << ' ' << zDir << std::endl;
    std::cout << "we are " << i << ' ' << j << std::endl;
    std::cout << x0 << ' ' << y0 << ' ' << z00 << std::endl;
    std::cout << x1 << ' ' << y0 << ' ' << z10 << std::endl;
    std::cout << x0 << ' ' << y1 << ' ' << z01 << std::endl;
    std::cout << x1 << ' ' << y1 << ' ' << z11 << std::endl;
    std::cout << z0 << ' ' << z1 << ' ' << z << std::endl;
    */

    return z;
}


