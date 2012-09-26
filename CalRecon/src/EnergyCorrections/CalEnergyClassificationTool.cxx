/**
 * @class CalEnergyClassificationTool
 *
 * @brief Implements a Gaudi Tool for selecting the best energy reconstruction method 
 *        from the choices available. This is driven off of a particular Tree/CalCluster
 *        association. 
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TreeClusterRelation.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"

#include "GlastSvc/GlastClassify/ITupleInterface.h"
#include "GlastSvc/GlastClassify/IClassifyTool.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "ICalEnergySelectionTool.h"
#include "CalUtil/IUBinterpolateTool.h"

#include <cmath>  // for M_PI, among others

// Define a concrete implementation of our "items"
// This is done here as a convenience in this module
template <class T> class ClassifyItem : public GlastClassify::Item 
{
public:
    ClassifyItem<T>(const std::string& name, const std::string& type, T* data) :
      m_pdata(data), m_name(name), m_type(type)   {}
    
    virtual ~ClassifyItem<T>() {}

    operator double()const {
        return (double)*m_pdata;
    }

    void*              getDataAddr() const {return m_pdata;}

    const std::string& getDataName() const {return m_name;}
    const std::string& getDataType() const {return m_type;}

// LSR 14-Jul-08 code for ntuple types

    void setDataValue(void* data) 
    {
        if (m_type == "UInt_t")
        {
            *m_pdata = *(reinterpret_cast<int*>(data));
        }
        else if (m_type == "ULong64_t")
        {
            *m_pdata = *(reinterpret_cast<unsigned long long*>(data));
        }
        else if (m_type == "Float_t")
        {
            *m_pdata = *(reinterpret_cast<float*>(data));
        }
        else if (m_type == "Double_t")
        {
            *m_pdata = *(reinterpret_cast<double*>(data));
        }
        else if (m_type == "UChar_t")
        {
            memset(m_pdata, ' ', 80);
            strcpy(reinterpret_cast<char*>(m_pdata), reinterpret_cast<char*>(data));
        }
        else if (m_type == "Char_t")
        {
            memset(m_pdata, ' ', 80);
            strcpy(reinterpret_cast<char*>(m_pdata), reinterpret_cast<char*>(data));
        }
        else
        {
            throw std::invalid_argument("ClassifyItem: attempting to set an unrecognized data type");
        }
    }

private:
    T*          m_pdata;
    std::string m_name;
    std::string m_type;
    void*       m_treePtr;
};


class CalEnergyClassificationTool : public AlgTool, virtual public ICalEnergySelectionTool
{
public:
    /// Standard Gaudi Tool interface constructor
    CalEnergyClassificationTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~CalEnergyClassificationTool() {}

    /// @brief This method initializes the tool, in particular setting up the classification tree
    ///        that will be used to select out the "best" energy for a particular tree/cluster pair
    StatusCode initialize();

    /// @brief This is the method called to run the classification on a particular tree/cluster pair
    ///        with the result then returned as a pointer to the best energy parameters
    const Event::CalCorToolResult* selectBestEnergy(Event::CalEventEnergy*      calEnergy,
                                                    Event::TreeClusterRelation* treeClusRel);


private:
    /// Internal methods
    bool   setTupleValues(Event::CalEventEnergy*      calEnergy,
                          Event::TreeClusterRelation* treeClusRel);

    /// This used to determine the active distance
    double CalEnergyClassificationTool::activeDistTower(Point pos, int &view) const;
    double CalEnergyClassificationTool::activeDistLAT(Point pos, int &view) const;

    /// Companion methods to above
    /// sign of a number
    static double sign(double x) { return x>0 ? 1.: -1. ;}

    /// turn a global coordinate (tower, ladder, wafer) roughly into a local one
    static double globalToLocal(double x, double pitch, int n) 
    {
        double xNorm = x / pitch + 0.5 * n;
        return sign(x) * (fmod(fabs(xNorm),1.0) - 0.5) * pitch ;
    }

    /// Define the root tuple file name and tuple/tree name here
    std::string                      m_tupleFileName;
    std::string                      m_tupleName;

    /// "Tuple" values used in this code
    IClassifyTool::VarNameToValueMap m_tupleMap;
    IClassifyTool::VarNameToValueMap m_outTupleMap;

    /// Local storage of variables
    std::map<std::string, float>     m_floatAddrMap;

    /// The output class
    char*                            m_outputString;

    /// Pointer to the ClassifyTool
    IClassifyTool*                   m_classifyTool;

    /// Minimum calorimeter energy
    double                           m_minCalEnergyRaw;

    /// xml analysis name 
    std::string                      m_analysisXmlFileName;

    /// Pointer to the Gaudi data provider service
    DataSvc*                         m_dataSvc;

    /// Unbiaser from CalUtil
    IUBinterpolateTool* m_ubInterpolateTool;

    /// some Geometry
    double m_towerPitch;
    int    m_xNum;
    int    m_yNum;
    //int    m_nLayers;
    //int    m_nCsI;
};

//static ToolFactory<CalEnergyClassificationTool> s_factory;
//const IToolFactory& CalEnergyClassificationToolFactory = s_factory;
DECLARE_TOOL_FACTORY(CalEnergyClassificationTool);

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

CalEnergyClassificationTool::CalEnergyClassificationTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ICalEnergySelectionTool>(this);

    // This allows for the outputting of the tuple that is used to feed the classificaiton tree to disk
    // If this variable is not null, then the tuple will be output to disk, if it is null the tuple
    // will be memory resident only. 
    // A typical format for the output name is something like: $(GLEAMDATAPATH)/CalEnergySelectionTuple.root
    declareProperty("TupleFileName",    m_tupleFileName       = "");

    // This allows one to define the name of the tree that would be output to disk (see above)
    declareProperty("TupleName",        m_tupleName           = "CalEnergySelection");

    // This defines the location and name of the xml file to be input for deciding which energy to use. 
    // A typical format for this variable is: $(CALRECONXMLPATH)/EnergyAnalysis_v2.xml
    declareProperty("AnalysisFileName", m_analysisXmlFileName = "$(CALRECONXMLPATH)/EnergyAnalysis_v2.xml");

    // Defines a minimum "raw" energy cut to run the classification
    declareProperty("MinCalEnergyRaw",  m_minCalEnergyRaw     = 10.);

    m_tupleMap.clear();
    m_outTupleMap.clear();
    m_floatAddrMap.clear();

    return;
}

StatusCode CalEnergyClassificationTool::initialize()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    // Take care of some needed geometry here...
    IGlastDetSvc* detSvc = 0;
    
    if (service("GlastDetSvc", detSvc, true).isFailure())
    {
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
        
    detSvc->getNumericConstByName("xNum",       &m_xNum);
    detSvc->getNumericConstByName("yNum",       &m_yNum);
    detSvc->getNumericConstByName("towerPitch", &m_towerPitch);


    // What follows below is the set up of the ntuple we use to feed to the classification analysis
    // While it appears to be a bit convoluted, it is straightforward:
    // 1) the map m_floatAddrMap maps the name of the variable to a float which will contain the actual variable
    // 2) "ClassifyItem" is a templated class that is stored in the tuple map, which maps the variable name to the "item"
    //    This structure allows us to store variables of multiple types in the same map, which is then fed to the 
    //    ClassifyTool to set up the processing of the classification analysis
    m_tupleMap["Cal1RawEnergySum"]     = new ClassifyItem<float>("Cal1RawEnergySum", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["Cal1RawEnergySum"]); // "Raw" energy sum cluster # 1?
    m_tupleMap["CalCsIRLn"]            = new ClassifyItem<float>("CalCsIRLn", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalCsIRLn"]);        // CsI radiation lengths this cluster
    m_tupleMap["CalTotRLn"]            = new ClassifyItem<float>("CalTotRLn", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalTotRLn"]);        // Total radiation lengths this tree/cluster
    m_tupleMap["CalEnergyRaw"]         = new ClassifyItem<float>("CalEnergyRaw", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalEnergyRaw"]);     // "Raw" energy this cluster
    m_tupleMap["CalLATEdge"]           = new ClassifyItem<float>("CalLATEdge", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalLATEdge"]);       // Closest distance to cal edge at top of cal
    m_tupleMap["CalLongRms"]           = new ClassifyItem<float>("CalLongRms", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalLongRms"]);       // Longitudinal rms of moments 
    m_tupleMap["CalTransRms"]          = new ClassifyItem<float>("CalTransRms", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalTransRms"]);      // transverse rms of moments
    m_tupleMap["CalTwrEdge"]           = new ClassifyItem<float>("CalLATEdge", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalLATEdge"]);       // Cal Tower edge
    m_tupleMap["CalTwrEdgeCntr"]       = new ClassifyItem<float>("CalTwrEdgeCntr", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalTwrEdgeCntr"]);   // Cal Tower edge center
    m_tupleMap["CalZDir"]              = new ClassifyItem<float>("CalLATEdge", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalLATEdge"]);       // cos(theta) for moments analysis
    m_tupleMap["CalCfpChiSq"]          = new ClassifyItem<float>("CalCfpChiSq", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalCfpChiSq"]);      // Profile fit chi square
    m_tupleMap["CalCfpEffRLn"]         = new ClassifyItem<float>("CalCfpEffRLn", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalCfpEffRLn"]);     // Profile fit effective radiation lengths
    m_tupleMap["CalCfpEnergy"]         = new ClassifyItem<float>("CalCfpEnergy", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalCfpEnergy"]);     // Profile fit energy
    m_tupleMap["CalCfpEnergyUB"]       = new ClassifyItem<float>("CalLATEdge", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalLATEdge"]);       // Profile fit energy
    m_tupleMap["CalCfpTkrRLn"]         = new ClassifyItem<float>("CalCfpTkrRLn", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["CalCfpTkrRLn"]);     // Profile fit radiation lengths in tracker
    m_tupleMap["EvtEnergyCorrUB"]      = new ClassifyItem<float>("EvtEnergyCorrUB", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["EvtEnergyCorrUB"]);  // Actually the energy from the parametric method
    m_tupleMap["Tkr1FirstLayer"]       = new ClassifyItem<float>("Tkr1FirstLayer", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["Tkr1FirstLayer"]);   // First layer of tree
    m_tupleMap["TkrTree1DirX"]         = new ClassifyItem<float>("TkrTree1DirX", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["TkrTree1DirX"]);     // direction cosine X for tree 1
    m_tupleMap["TkrTree1DirY"]         = new ClassifyItem<float>("TkrTree1DirY", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["TkrTree1DirY"]);     // direction cosine Y for tree 1
    m_tupleMap["TkrTree1DirZ"]         = new ClassifyItem<float>("TkrTree1DirZ", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["TkrTree1DirZ"]);     // direction cosine Z for tree 1


    m_tupleMap["McEnergy"]             = new ClassifyItem<float>("McEnergy", 
                                                                 "Float_t", 
                                                                 &m_floatAddrMap["McEnergy"]);         // McEnergy but should not be here

    // The end is near! 
    // Recover the ClassifyTool from the land of Gaudi
    if ((sc = toolSvc()->retrieveTool("ClassifyTool", "ClassifyTool", m_classifyTool)).isFailure())
    {
        throw GaudiException("Service [ClassifyTool] not found", name(), sc);
    }

    // And intitialize it
    m_classifyTool->setUpClassification(m_tupleMap, m_analysisXmlFileName, m_tupleName, m_tupleFileName);

    // The output of the classification is handled differently from the intput, in that the "Item"'s for 
    // the output are set up in the above initialization. To keep local copies we simply look up those
    // Item's and keep track of them in our own internal map
    GlastClassify::Item* energyMethod = 0;

    if (m_classifyTool->getVariable("CTBEnergyMethod", energyMethod))
    {
        m_outTupleMap["CTBEnergyMethod"] = energyMethod;
    }

    // Finally get a hold of the unbiasing interpolation code
    // Get the Tool
    m_ubInterpolateTool = 0;
    // Recover the ClassifyTool from the land of Gaudi
    if ((sc = toolSvc()->retrieveTool("UBinterpolateTool", "UBinterpolateTool", m_ubInterpolateTool)).isFailure())
    {
        throw GaudiException("Service [UBinterpolateTool] not found", name(), sc);
    }
    m_ubInterpolateTool->addBiasMap("Prof","$(CALUTILXMLPATH)/BiasMapCalCfpEnergy.txt"); 
    //m_ubInterpolate = new UBinterpolate("$(CALRECONXMLPATH)/BiasMapCalCfpEnergy.txt");
    
    return sc;
}

const Event::CalCorToolResult* CalEnergyClassificationTool::selectBestEnergy(Event::CalEventEnergy*      calEnergy,
                                                                             Event::TreeClusterRelation* treeClusRel)
{
    // Purpose and Method: Determine the "best" energy from the available choices
    // Inputs:  Tree/Cluster relation and the list of available reconstructed energies
    // Outputs: The object containing the parameters of the "best" energy 
    // Dependencies: None
    // Restrictions and Caveats:  None.
    //Always believe deep down in the success of this venture we are about to embark on
    StatusCode sc = StatusCode::SUCCESS;

    // Define a pointer to the "corrected" energy result
    const Event::CalCorToolResult* calResult = 0;

    // Check that we have some objects in place...
    if (calEnergy && treeClusRel)
    {
        // Set up the tuple to use here
        if (setTupleValues(calEnergy, treeClusRel))
        {
            // Call the classification routine
            m_classifyTool->runClassification();

            // Get back the result
            std::string energyMethod((char*)m_outTupleMap["CTBEnergyMethod"]->getDataAddr());

            if (energyMethod == "Profile") calResult = calEnergy->findLast("CalFullProfileTool");
            else                           calResult = calEnergy->findLast("CalValsCorrTool");
        }
    
        // The default value, in case nothing happens, is to return the "raw" energy
        if (!calResult) calResult = calEnergy->findLast("CalRawEnergyTool");
    }

    return calResult;
}
    
bool CalEnergyClassificationTool::setTupleValues(Event::CalEventEnergy*      calEnergy,
                                                 Event::TreeClusterRelation* treeClusRel)
{
    // Set the "tuple" values that are used in this objects classification tree here
    // First, dereference the tree and cluster
    Event::TkrTree*    tree    = treeClusRel->getTree();
    Event::CalCluster* cluster = treeClusRel->getCluster();

    // No tree, no cluster, no work
    if (!cluster) return false;

    // useful stuff
    Vector treeDir(0., 0., -1.);
    float  tkrDirZ = treeDir.z();

    // Recover the energy corrections here
    const Event::CalCorToolResult* calValsCorr = calEnergy->findLast("CalValsCorrTool");
    const Event::CalCorToolResult* calCfpCorr  = calEnergy->findLast("CalFullProfileTool");

    // if no energy corrections then we should return false
    if (!calValsCorr) return false;

    // CalEnergyRaw 
    Point  calPos       = cluster->getPosition();
    double calEnergyRaw = cluster->getXtalsParams().getXtalRawEneSum();
    int    tmpView      = 0;

    // ugliness to calculate the two edge variables we need
    float twrEdgeTop = 0.;
    float latEdge    = 0.;
    
    if (tree)
    {
        // start with hardwired number for cal z top...
        double calZTop = -48.12 - 9.5;

        // Tree parameters
        Point  treePos = tree->getAxisParams()->getEventPosition();
        
        treeDir = -tree->getAxisParams()->getEventAxis();
        tkrDirZ = treeDir.z();

        // Distance to top of Cal
        double arcLenToTop = (treePos.z() - calZTop) / treeDir.z();

        // Position at the top of the Cal
        Point calTopPos = treePos - arcLenToTop * treeDir;

        // Get the activedistance
        int tmpView = 0;

        // Ok, first variable we need is the active distance to nearest edge
        twrEdgeTop = activeDistTower(calTopPos, tmpView);

        // I think the next is meant to be the distance to the edge of the LAT
        latEdge = activeDistLAT(calTopPos, tmpView);
    }

    calEnergyRaw = std::max(m_minCalEnergyRaw, calEnergyRaw);

    m_floatAddrMap["CalEnergyRaw"]     = calEnergyRaw;

    m_floatAddrMap["Cal1RawEnergySum"] = cluster->getXtalsParams().getXtalRawEneSum(); 

    m_floatAddrMap["CalZDir"]          = cluster->getMomParams().getAxis().z();
    m_floatAddrMap["CalLongRms"]       = cluster->getMomParams().getLongRms();
    m_floatAddrMap["CalTransRms"]      = cluster->getMomParams().getTransRms();

    m_floatAddrMap["CalLATEdge"]       = latEdge;
    m_floatAddrMap["CalTwrEdge"]       = twrEdgeTop;
    m_floatAddrMap["CalTwrEdgeCntr"]   = activeDistTower(calPos, tmpView);


    if (calCfpCorr)
    {
        float calCfpEnergy = calCfpCorr->getParams().getEnergy();

        m_floatAddrMap["CalCfpEnergy"]   = calCfpEnergy;
        m_floatAddrMap["CalCfpChiSq"]    = calCfpCorr->find("totchisq")->second;
        m_floatAddrMap["CalCfpEffRLn"]   = calCfpCorr->find("cal_eff_RLn")->second;
        m_floatAddrMap["CalCfpTkrRLn"]   = calCfpCorr->find("tkr_RLn")->second;

        // need to do the event unbiasing here
        const float bias = m_ubInterpolateTool->interpolate("Prof", log10(calCfpEnergy), tkrDirZ);

        float calCfpEnergyUB = bias == 0 ? -1 : calCfpEnergy / bias;
        
        m_floatAddrMap["CalCfpEnergyUB"] = calCfpEnergyUB;
    }
    else
    {
        m_floatAddrMap["CalCfpEnergy"]   = 0.;
        m_floatAddrMap["CalCfpChiSq"]    = 0.;
        m_floatAddrMap["CalCfpEffRLn"]   = 0.;
        m_floatAddrMap["CalCfpEnergyUB"] = 0.;
        m_floatAddrMap["CalCfpTkrRLn"]   = 0.;
    }

    // THIS IS NOT RIGHT - need to recover the unbiased event energy here... but don't know how 
    // to do here, leave for the experts to determine how to include
    m_floatAddrMap["EvtEnergyCorrUB"] = calValsCorr->find("CorrectedEnergy")->second; 
    m_floatAddrMap["CalCsIRLn"]       = calValsCorr->find("CsIRLn")->second;
    m_floatAddrMap["CalTotRLn"]       = calValsCorr->find("CALRLn")->second;

    // Recover tree vars if there is a tree
    if (tree)
    {
        m_floatAddrMap["Tkr1FirstLayer"] = tree->getHeadNode()->getTreeStartLayer();
        m_floatAddrMap["TkrTree1DirX"]   = treeDir.x();
        m_floatAddrMap["TkrTree1DirY"]   = treeDir.y();
    }
    else
    {
        m_floatAddrMap["Tkr1FirstLayer"] = 0.;
        m_floatAddrMap["TkrTree1DirX"]   = 0.;
        m_floatAddrMap["TkrTree1DirY"]   = 0.;
    }

    m_floatAddrMap["TkrTree1DirZ"]         = tkrDirZ;

    m_floatAddrMap["McEnergy"]             = 0.;

    return true;
}

double CalEnergyClassificationTool::activeDistTower(Point pos, int &view) const
{
    double edge = 0.;
    double x = pos.x();
    double y = pos.y();
    double x_twr = globalToLocal(x, m_towerPitch, m_xNum);
    double y_twr = globalToLocal(y, m_towerPitch, m_yNum);

    if( fabs(x_twr) > fabs(y_twr) ) {
        edge = m_towerPitch/2. - fabs(x_twr);
        view = 0; 
    }
    else {
        edge = m_towerPitch/2. - fabs(y_twr);
        view = 1;
    }
    return edge;
}

double CalEnergyClassificationTool::activeDistLAT(Point pos, int &view) const
{
    double edge    = 0.;
    double latEdge = (m_xNum * m_towerPitch) / 2.;
    double x_lat   = pos.x();
    double y_lat   = pos.y();

    if( fabs(x_lat) > fabs(y_lat) ) {
        edge = latEdge - fabs(x_lat);
        view = 0; 
    }
    else {
        edge = latEdge - fabs(y_lat);
        view = 1;
    }
    return edge;
}
