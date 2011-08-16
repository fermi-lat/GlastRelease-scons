
// Tool and Gaudi related stuff
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalNBCClassParams.h"

#include <CalRecon/ICalReconSvc.h>
#include <CalRecon/ICalClassifyTool.h>

// To parse the xml file with PDFs
#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "facilities/Util.h"
#include "facilities/commonUtilities.h"

#include <algorithm>
#include <iostream>


class NBCEnergyBin
{
public:
    NBCEnergyBin(double logEmin, double logEmax):
      m_logEmin(logEmin), m_logEmax(logEmax) {;}
      ~NBCEnergyBin() {};

      inline double getLoEdge() const {return m_logEmin;}
      inline double getHiEdge() const {return m_logEmax;}
      inline double getCenter() const {return 0.5*(getLoEdge() + getHiEdge());}

      friend std::ostream& operator<<(std::ostream& output, const NBCEnergyBin& bin);

private:
    double m_logEmin;
    double m_logEmax;
};

std::ostream& operator<<(std::ostream& output, const NBCEnergyBin& bin)
{
    output << "log(E [MeV]) = (" << bin.getLoEdge() << "--" <<
        bin.getHiEdge() << "), center @ " << bin.getCenter();
    return output;
}



class NBCEnergyBinning
{
public:
    NBCEnergyBinning() {;}
    ~NBCEnergyBinning() {};

    inline int getNumBins()              const {return m_energyBins.size();}
    inline NBCEnergyBin getBin(int bin)  const {return m_energyBins[bin];}
    inline NBCEnergyBin getFirstBin()    const {return getBin(0);}
    inline NBCEnergyBin getLastBin()     const {return getBin(getNumBins() - 1);}
    inline double getBinLoEdge(int bin)  const {return getBin(bin).getLoEdge();}
    inline double getMinimum()           const {return getFirstBin().getLoEdge();}
    inline double getBinHiEdge(int bin)  const {return getBin(bin).getHiEdge();}
    inline double getMaximum()           const {return getLastBin().getHiEdge();}
    inline double getBinCenter(int bin)  const {return getBin(bin).getCenter();}

    void addBin(double logEmin, double logEmax);
    int getBinIndex(double energy);

    friend std::ostream& operator<<(std::ostream& output,
        const NBCEnergyBinning& binning);

private:
    std::vector <NBCEnergyBin> m_energyBins;
};

void NBCEnergyBinning::addBin(double logEmin, double logEmax)
{
    m_energyBins.push_back(NBCEnergyBin(logEmin, logEmax));
}

int NBCEnergyBinning::getBinIndex(double energy)
{
    if (energy <= 0) return 0;
    double logE = log10(energy);
    if (logE <= getMinimum()) return 0;
    for (int i = 0; i < getNumBins(); i++)
    {
        if (logE < getBin(i).getHiEdge()) return i;
    }
    return getNumBins() - 1;
}

std::ostream& operator<<(std::ostream& output, const NBCEnergyBinning& binning)
{
    for (int i = 0; i < binning.getNumBins(); i++)
    {
        output << "Bin " << i << ": " << binning.getBin(i) << std::endl;
    }
    return output;
}


/*
Tiny class to store histograms with the PDF
for one variable for one energy bin
*/

class PdfHisto
{
public:
    PdfHisto() {;}
    ~PdfHisto() {};

    inline int getNumBins()             const {return m_probs.size();}
    inline double getBinLoEdge(int bin) const {return m_loEdges[bin];}
    inline double getBinHiEdge(int bin) const {return m_hiEdges[bin];}
    inline double getBinProb(int bin)   const {return m_probs[bin];}

    void addBin(double xmin, double xmax, double prob);
    double getPdValue(double varValue);

    friend std::ostream& operator<<(std::ostream& output, const PdfHisto& hist);

private:
    std::vector <double> m_loEdges;
    std::vector <double> m_hiEdges;
    std::vector <double> m_probs;
};

void PdfHisto::addBin(double xmin, double xmax, double prob)
{
    m_loEdges.push_back(xmin);
    m_hiEdges.push_back(xmax);
    m_probs.push_back(prob);
}

// Sorting?
// Uderflow/overlfow?
double PdfHisto::getPdValue(double varValue)
{   
    std::vector<double>::iterator it;
    int index=0;
    for (it=m_hiEdges.begin() ; it != m_hiEdges.end(); it++ )
    {
        if (varValue<(*it))
        {
            return m_probs[index];
        }
        index++;
    }
    return 0;
}

std::ostream& operator<<(std::ostream& output, const PdfHisto& hist)
{
    for (int i = 0; i < hist.getNumBins(); i++)
    {
        output << "Bin " << i << ": (" << hist.getBinLoEdge(i) << "--" << hist.getBinHiEdge(i) 
            << "), prob = " << hist.getBinProb(i) << std::endl;
    }
    return output;
}



/*
Collection of PdfHisto objects
one PdfHisto per energy bin
one collection for each variable for each topology
*/

class PdfHistoCollection
{
public:
    PdfHistoCollection(std::string topo, std::string varName, NBCEnergyBinning *energyBinning):
      m_topology(topo), m_varName(varName), m_energyBinning(energyBinning) {;}
      ~PdfHistoCollection() {};

      inline std::string getTopology()           {return m_topology;}
      inline std::string getVarName()            {return m_varName;}
      inline void addPdfHisto(PdfHisto hist)     {m_pdfHistos.push_back(hist);}

      PdfHisto getPdfHisto(double energy);
      double getPdValue(double energy, double varValue);

private:
    std::string m_topology;
    std::string m_varName;
    NBCEnergyBinning *m_energyBinning;
    std::vector <PdfHisto> m_pdfHistos;
};

PdfHisto PdfHistoCollection::getPdfHisto(double energy)
{
    return m_pdfHistos[m_energyBinning->getBinIndex(energy)];
}

double PdfHistoCollection::getPdValue(double energy, double varValue)
{
    return getPdfHisto(energy).getPdValue(varValue);
}

/*
Definition of the tool for the NBC classifier
*/

class CalClusterNBClassifyTool : public AlgTool, virtual public ICalClassifyTool
{
public :

    /// Standard Gaudi Tool interface constructor
    CalClusterNBClassifyTool(const std::string& type,
        const std::string& name,
        const IInterface* parent );

    virtual ~CalClusterNBClassifyTool() {};

    /// @brief Intialization of the tool
    virtual StatusCode initialize() ;

    /// @brief Default cluster finding framework
    virtual StatusCode classifyClusters(Event::CalClusterCol* calClusterCol);

    virtual StatusCode classifyCluster(Event::CalCluster* calCluster);

private:
    //! Service for basic Cal info
    ICalReconSvc*      m_calReconSvc;

    //! Event Service member directly useable by concrete classes.
    IDataProviderSvc*  m_dataSvc;

    // Other methods
    StatusCode getPDFsFromXml();
    double getPdValue(std::string topology, std::string varName,
        double energy, double varValue);
    double getVariableValue(std::string varName, Event::CalCluster* calCluster);
    void addPdfHistCollection(PdfHistoCollection pdfHistCol);
    PdfHistoCollection getHistoCollection(std::string topology,
        std::string varName);

    // Members
    std::string m_xmlPDFFileName;
    std::map <std::pair <std::string, std::string>, PdfHistoCollection> m_varPDFsMap;
    NBCEnergyBinning* m_energyBinning;    
};

DECLARE_TOOL_FACTORY(CalClusterNBClassifyTool) ;

CalClusterNBClassifyTool::CalClusterNBClassifyTool(const std::string & type, 
                                                   const std::string & name, 
                                                   const IInterface* parent)
                                                   : AlgTool(type,name,parent)
{ 
    declareInterface<ICalClassifyTool>(this) ;

    
    declareProperty ("m_xmlPDFFileName", m_xmlPDFFileName = "$(CALRECONXMLPATH)/xml_cluclass_NIPsWoTwoVars_CalTwrEdge.xml" );

    return;
}

void CalClusterNBClassifyTool::addPdfHistCollection(PdfHistoCollection pdfHistCol)
{
    std::pair <std::string, std::string> key(pdfHistCol.getTopology(), pdfHistCol.getVarName());
    m_varPDFsMap.insert(std::make_pair(key, pdfHistCol));
}

StatusCode CalClusterNBClassifyTool::initialize()
{
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc = StatusCode::SUCCESS;

    if ((sc = service("CalReconSvc",m_calReconSvc,true)).isFailure())
    {
        throw GaudiException("Service [CalReconSvc] not found", name(), sc);
    }

    if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    // check if jobOptions were passed
    log<<MSG::INFO<<"m_xmlPDFFileName\t\tset to "<<m_xmlPDFFileName<<endreq;

    // Cluster filling utility
    // TBD -- Need to update cluster with a classification info
    // m_clusterInfo = new MomentsClusterInfo(m_calReconSvc);

    // Read the xml file with PDFs
    if ( sc = getPDFsFromXml().isFailure() )
        log<<MSG::ERROR<<"Could not read xml file"<<endreq;

    return StatusCode::SUCCESS ;
}

StatusCode CalClusterNBClassifyTool::classifyClusters(Event::CalClusterCol* calClusterCol)
{
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc = StatusCode::SUCCESS;
    //Purpose and method:
    //
    //   This function performs the classification of the calorimeter clusters
    //   The main actions are:
    //      - classify the cluster
    // TDS output: updated CalClustersCol
    // -- TBD --

    // --------------------------------------------

    log << MSG::DEBUG << "Calling Naive Bayes classifier" << endreq;
    Event::CalClusterCol::const_iterator cluster ;
    for (cluster = calClusterCol->begin(); cluster != calClusterCol->end(); cluster++)
    {
        // Current PDFs were derived with precut at NumXtals>3
        if( (*cluster)->getNumXtals() > 3 ) {
            classifyCluster(*cluster);
        }
    }   

    return StatusCode::SUCCESS ;
}

StatusCode CalClusterNBClassifyTool::classifyCluster(Event::CalCluster* calCluster)
{
    MsgStream log(msgSvc(),name()) ;
    std::string className, varName;
    double varValue, pdValue;
    double energy = calCluster->getMomParams().getEnergy();
    std::map<std::pair <std::string, std::string>, PdfHistoCollection>::iterator iter;

    // Create the CalNBCClassParams object to store the class probabilities.
    Event::CalNBCClassParams nbcClassParams;
    for (iter = m_varPDFsMap.begin(); iter != m_varPDFsMap.end(); iter++)
    {
        className = (*iter).first.first;
        varName = (*iter).first.second;
        varValue = getVariableValue(varName, calCluster);
        pdValue = getPdValue(className, varName, energy, varValue);
        nbcClassParams.multiply(className, pdValue);
    }
    nbcClassParams.normalize();

    // And finally set the corresponding cluster class member.
    calCluster->setClassParams(nbcClassParams);
    calCluster->setStatusBit(Event::CalCluster::CLASSIFIED);

    return StatusCode::SUCCESS;
}

double CalClusterNBClassifyTool::getVariableValue(std::string varName,
                                                  Event::CalCluster* calCluster)
{
    if (varName == "CalTransRms") {
        return calCluster->getMomParams().getTransRms();
    }
    else if (varName == "CalLongRmsAsym") {
        return calCluster->getMomParams().getLongRmsAsym();
    }
    else if (varName == "CalCoreEneFrac") {
        return calCluster->getMomParams().getCoreEnergyFrac();
    }
    else if (varName == "CalMomNumCoreXtals") {
        return calCluster->getMomParams().getNumCoreXtals();
    }
    else if (varName == "XtalEneRms") {
      double sumOfWeights = calCluster->getMomParams().getEnergy();
      if (sumOfWeights > 0)
        return (calCluster->getXtalsParams().getXtalEneRms())/sumOfWeights;
      else
	return -1;
    }
    else if (varName == "CalXtalEneSkewness") {
        return calCluster->getXtalsParams().getXtalEneSkewness();
    }
    else if (varName == "CalMstAveEdgeLen") {
        return calCluster->getMSTreeParams().getMeanEdgeLength();
    }
    else if (varName == "dEdxperLength") {
      double MomdEdxAve  = calCluster->getMomParams().getdEdxAverage();
      double MomFullLenght = calCluster->getMomParams().getFullLength();
      return log10(MomdEdxAve/MomFullLenght);     
    }
    else if (varName == "CalLyr0Ratio") {
      double CALEnergyRaw  = calCluster->getMomParams().getEnergy();
      double ELayer0 = calCluster->at(0).getEnergy();
      if (CALEnergyRaw > 0.0) 
	return ELayer0/CALEnergyRaw;
      else
	return -1;
    }
    else if (varName == "MomentRatio") {
        double transRms = calCluster->getMomParams().getTransRms();
        double longRms = calCluster->getMomParams().getLongRms();
        if ( transRms > 0 && longRms > 0 ) {
            return log10(longRms/transRms);
        }
        else {
            return -1.;
        }
    }
    else if (varName == "NumXtals") {
        double sumOfWeights = calCluster->getMomParams().getEnergy();
        if (sumOfWeights > 0)
            return (calCluster->getNumXtals())/log10(sumOfWeights);
        else
            return -1.0;
    }
    else {
        std::cout << "Unknown variable! NBC will not work!" << varName << ". Abort." << std::endl;
        return -1.0;
    }
}

StatusCode CalClusterNBClassifyTool::getPDFsFromXml()
{
    XERCES_CPP_NAMESPACE_USE
        MsgStream log(msgSvc(),name());
    facilities::commonUtilities::setupEnvironment();

    facilities::Util::expandEnvVar(&m_xmlPDFFileName);

    xmlBase::XmlParser* parser = new xmlBase::XmlParser(true);

    DOMDocument* doc = 0;
    try {
        doc = parser->parse(m_xmlPDFFileName.c_str());
    }
    catch (xmlBase::ParseException ex) {
        std::cout << "caught exception with message " << std::endl;
        std::cout << ex.getMsg() << std::endl;
        delete parser;
        return StatusCode::FAILURE;
    }

    if (doc != 0)   // successful
    {
        log << MSG::DEBUG << "Document successfully parsed"  << endreq;

        // look up some attributes
        DOMElement* docElt = doc->getDocumentElement();

        // Some variables used everywhere
        int    intVal;

        // ----------------------------------
        // Read the first tag with energy bins
        // ----------------------------------
        m_energyBinning = new NBCEnergyBinning();

        DOMElement* attElt = xmlBase::Dom::findFirstChildByName(docElt, "EnergyBins");       

        std::vector<DOMElement*> childrenEnergyBins;
        std::vector<DOMElement*>::iterator itemElt;
        // Get a vector of DOMElements with Tag Name Energy, corresponding to all energy bins
        xmlBase::Dom::getChildrenByTagName(attElt, "Energy",childrenEnergyBins,true);
        for (itemElt=childrenEnergyBins.begin() ; itemElt != childrenEnergyBins.end(); itemElt++ )
        {
            double emin = xmlBase::Dom::getDoubleAttribute(*itemElt, "Emin");
            double emax = xmlBase::Dom::getDoubleAttribute(*itemElt, "Emax");
            m_energyBinning->addBin(emin, emax);
        }

        log<<MSG::DEBUG << "Energy binning for the NBC classifier read from xml" <<
            std::endl << *m_energyBinning << endreq;

        // ----------------------------------
        // Get the vector of topologies DOMElement
        // ----------------------------------       
        double xmin, xmax, prob;
        std::vector<DOMElement*> topoVector;
        xmlBase::Dom::getChildrenByTagName(docElt, "Topology", topoVector, true);
        log<<MSG::DEBUG << topoVector.size() << 
            " cluster topologies found in the xml file." << endreq; 

        // loop on topologies
        std::vector<DOMElement*>::iterator topoIter;                
        for (topoIter=topoVector.begin() ; topoIter != topoVector.end(); topoIter++ ) 
        {            
            // Get the name of the topology
            std::string topology = xmlBase::Dom::getAttribute(*topoIter, "name");

            // Get the vector of variables DOMElement
            std::vector<DOMElement*> varVector;
            xmlBase::Dom::getChildrenByTagName(*topoIter, "Variable", varVector, true);
            log<<MSG::DEBUG << varVector.size() << " variables defined for topology '"
                << topology <<"'." << endreq;

            // loop on variables
            std::vector<DOMElement*>::iterator varIter;                
            for (varIter=varVector.begin() ; varIter != varVector.end(); varIter++ ) 
            {            
                // Get the name of the variable
                std::string varName = xmlBase::Dom::getAttribute(*varIter, "name");

                // Create the PdfHisoCollection that we'll fill looping on energy bins
                log << MSG::DEBUG << "Instantiating pdf collection for topology '" <<
                    topology << "' and variable '" << varName << "'" << endreq;
                PdfHistoCollection pdfHistCol(topology, varName, m_energyBinning);

                // Get the vector of Energy bins DOMElement
                std::vector<DOMElement*> eneVector;
                xmlBase::Dom::getChildrenByTagName(*varIter, "Energy", eneVector, true);

                // loop on Energy bins for the variable
                std::vector<DOMElement*>::iterator eneIter;                 
                for (eneIter=eneVector.begin() ; eneIter != eneVector.end(); eneIter++ ) 
                {
                    PdfHisto hist;            
                    // Get the value of the energy bin for the variable
                    intVal = xmlBase::Dom::getIntAttribute(*eneIter, "bin");

                    // Get the vector of BinValues bins DOMElement -- Data !
                    std::vector<DOMElement*> dataVector;
                    xmlBase::Dom::getChildrenByTagName(*eneIter, "BinValues", dataVector, true);

                    // loop on Data bins BinValues for the Energy bin
                    std::vector<DOMElement*>::iterator dataIter;                 
                    for (dataIter=dataVector.begin() ; dataIter != dataVector.end(); dataIter++ )
                    {
                        xmin = xmlBase::Dom::getDoubleAttribute(*dataIter, "xmin");
                        prob = xmlBase::Dom::getDoubleAttribute(*dataIter, "pdv");
                        xmax = xmlBase::Dom::getDoubleAttribute(*dataIter, "xmax");
                        hist.addBin(xmin, xmax, prob);
                    } // end of loop on Data
                    pdfHistCol.addPdfHisto(hist);
                    //log << MSG::DEBUG <<hist<<endreq;
                }// end of loop on  Energy 
                addPdfHistCollection(pdfHistCol);
            }// end of loop on  Variable
        }// end of loop on  Topology
    } // end of if parsing ok

    delete parser;
    return StatusCode::SUCCESS ;
}

PdfHistoCollection CalClusterNBClassifyTool::getHistoCollection(std::string topology,
                                                                std::string varName)
{
    return m_varPDFsMap.find(std::make_pair(topology, varName))->second;
}

double CalClusterNBClassifyTool::getPdValue(std::string topology, std::string varName,
                                            double energy, double varValue)
{
    return getHistoCollection(topology, varName).getPdValue(energy, varValue);
}
