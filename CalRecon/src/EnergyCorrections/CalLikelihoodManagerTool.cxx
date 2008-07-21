
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "facilities/Util.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "CalLikelihoodManagerTool.h"
#include "CalLikelihoodPDFFunctions.h"
#include <string>

/******************************************************************************/
/***********************  CalLikelihoodManagerTool ****************************/
/******************************************************************************/
#include <GaudiKernel/DeclareFactoryEntries.h>
DECLARE_TOOL_FACTORY(CalLikelihoodManagerTool) ;

CalLikelihoodManagerTool::CalLikelihoodManagerTool(const std::string& type,
                                                   const std::string& name, 
                                                   const IInterface* parent)
  : AlgTool(type,name,parent), PDFAnswer()
{
    // declare base interface for all consecutive concrete classes
    std::vector<std::string> parameters;
		parameters.push_back("$(CALRECONXMLPATH)/HighEnergyLowIncidence.data");
		parameters.push_back("$(CALRECONXMLPATH)/LowEnergy.data");
		parameters.push_back("$(CALRECONXMLPATH)/HighEnergyHighIncidence.data");

    
    declareInterface<ICalEnergyCorr>(this);
    declareProperty("parameterFiles", m_parameterFiles= parameters);
}
StatusCode CalLikelihoodManagerTool::initialize(void)
{
  // This function does following initialization actions:
  //    - extracts geometry constants from xml file using GlastDetSvc
  MsgStream log(msgSvc(), "CalLikelihoodManagerTool::initialise");
  StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing"<<endreq;

  //Locate and store a pointer to the data service
  // which allows access to the TDS
  if ((sc = service("EventDataSvc", m_dataSvc)).isFailure())
  {
      throw GaudiException("Service [EventDataSvc] not found", name(), sc);
  }

  if ((sc = service("GlastDetSvc", m_detSvc, true)).isFailure())
  { 
    throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
  }

  //Get geometry info
  m_numX = 4;
  m_numY = 4;
  m_flight_geom = true;
  if(!m_detSvc->getNumericConstByName("xNum", &m_numX))
    {
      log<< MSG::ERROR << "constant " << "xNum" << " not defined" <<endreq;
      throw GaudiException("Bad constant name ", name(),StatusCode::FAILURE);
    }
  if(!m_detSvc->getNumericConstByName("yNum", &m_numY))
    {
      log<< MSG::ERROR << "constant " << "yNum" << " not defined" <<endreq;
      throw GaudiException("Bad constant name ", name(),StatusCode::FAILURE);
    }
  if(m_numX==4 && m_numY==1) m_flight_geom = false;

  log << MSG::INFO << "m_numX = " << m_numX << endreq;  
  log << MSG::INFO << "m_numY = " << m_numY << endreq;  
  if(m_flight_geom)
    log << MSG::INFO << "flight mode" << endreq;  
  else
    log << MSG::INFO << "calibration unit mode" << endreq;  

  m_minTrialEnergy= 1.e30;
  m_maxTrialEnergy= -1.e30;
  m_minTkr1ZDir= 2.;

  int nFiles= 0;
  const std::vector< std::string >& parametersFiles = m_parameterFiles ;
  std::vector< std::string >::const_iterator file;
  
  for(file= parametersFiles.begin(); file!=parametersFiles.end(); ++file) {
    PDFCutArray *tool= new PDFCutArray((*file), this, log);
    if( !(tool[0])  ) delete tool;
    else {
      ++nFiles;
      PDFVect::push_back(tool);
      PDFGrid *axes= tool->getGrid();
      double tmp= axes->getBinCenter(0, 0);
      if( tmp<m_minTrialEnergy ) m_minTrialEnergy= tmp;
      tmp= axes->getBinCenter(0, -1);
      if( tmp>m_maxTrialEnergy ) m_maxTrialEnergy= tmp;
      tmp= axes->getBinCenter(1, 0);
      if( tmp<m_minTkr1ZDir ) m_minTkr1ZDir= tmp;
      continue;
    }

    log<<MSG::ERROR<<"PDFCutArray could not be created."<<endreq
                   <<"At PDFCutArray: "<<nFiles<<endreq;
    clear();
    break;
  }    

  if( !nFiles ){
    log<<MSG::ERROR<<"CalLikelihoodManagerTool: "
                     "No PDFCutArray Created"<<endreq;
    clear();
  }
  return operator!()?StatusCode::FAILURE:StatusCode::SUCCESS;
}

StatusCode CalLikelihoodManagerTool::finalize(void)
{
  StatusCode sc = StatusCode::SUCCESS;
  clear();
  return sc;
}


void CalLikelihoodManagerTool::getParameters(
                                std::map<double*,std::string> &param,
                                MsgStream &log) const {
  
  for(std::map<double*,std::string>::iterator it=param.begin(); it!=param.end();
      it++)
  {
      if(!m_detSvc->getNumericConstByName((*it).second, (*it).first)) 
      {
          log<<MSG::ERROR<<"constant "<<(*it).second<<" not defined"<<endreq;
          throw GaudiException("Bad constant name ", name(),
                                StatusCode::FAILURE);
      }
  }
}
const Event::TkrDigiCol *CalLikelihoodManagerTool::getTkrDigiCol(void) const
{ 
  SmartDataPtr<Event::TkrDigiCol> digiCol(m_dataSvc, 
                                          EventModel::Digi::TkrDigiCol );
  return digiCol;
}

Event::CalCorToolResult * CalLikelihoodManagerTool::doEnergyCorr
 ( Event::CalClusterCol * clusters, Event::TkrVertex* vertex ){
                             
  MsgStream log(msgSvc(), "CalLikelihoodManagerTool.doEnergyCorr");
  if (vertex == 0)
  {
    log << MSG::DEBUG <<"No TKR Reconstruction"<< endreq;
    return 0;
  }

  if (clusters->empty())
  {
      log << MSG::DEBUG << "Ending doEnergyCorr: No Cluster" 
          << endreq;
      return 0 ;
  }
  Event::CalCluster * cluster = clusters->front() ;

  double tkr1Zdir= -vertex->getDirection()[2];
  double eMin= cluster->getCalParams().getEnergy();
  double eMax= eMin*5<100.?100.:eMin*5;

  if( getTkrPlane(vertex)<0 
      || eMax<m_minTrialEnergy 
      || eMin>m_maxTrialEnergy
      || tkr1Zdir<m_minTkr1ZDir) {
    log<<MSG::DEBUG<<"Outside range with CalEnergyRaw= "<<eMin<<" MeV, "
                     "Tkr1ZDir= "<<-tkr1Zdir
                   <<", Vertex #"<<getTkrPlane(vertex)<<endreq;
    return 0;
  }
  if( evalEnergy(eMin, eMax, cluster, vertex) || evalError(vertex) ) {
    log<<MSG::DEBUG<<"No PDF for this event"<<endreq;
    return 0;
  }

  Event::CalCorToolResult *result= new Event::CalCorToolResult();
  result->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
  result->setCorrectionName(type());
  result->setChiSquare(1.);
  Event::CalParams params= cluster->getCalParams();
  params.setEnergy(getEnergy());
  params.setEnergyErr(getError());
  result[0]["Probability"]= getProbability();
  result[0]["ErrorLow"]= getError(0);
  result[0]["ErrorHigh"]= getError(1);
  result[0]["AnswerID"]= getAnswerID()*10+getAnswer()->getAnswerID();

  result[0]["TkrSumHits"]= ((PDFLikelihood*) ((PDFCutArray*) getAnswer())
                            ->getLikelihoodFcn())
                            ->tkrSumHits();
  // next line means that the second PDFCutArray should be a "LowEnergy" type
  result[0]["GeometricCut"]= ((PDFLowEnergyCuts*) ((PDFCutArray*) at(1))
                              ->getCutsFcn())
                              ->geometricCut(cluster, vertex);
  result->setParams(params);

  log<<MSG::DEBUG<<"Energy "<<params.getEnergy()
                 <<" MeV, Error "<<params.getEnergyErr()
                 <<" Range "<<result[0]["AnswerID"]<<endreq;
  return result;
}

/******************************************************************************/
/***********************  PDFCutArray *****************************************/
/******************************************************************************/
PDFCutArray::PDFCutArray(std::string file,
                         const CalLikelihoodManagerTool* tool, 
                         MsgStream &log)
    : PDFAnswer(), m_Narray(0), m_Grid(0),
      m_LikelihoodFcn(new PDFLikelihood(tool, log)), 
      m_CutsFcn(0), m_Cuts(0){
  facilities::Util::expandEnvVar(&file);
  std::ifstream dataFile(file.data());
  if( !dataFile ){
    log<<MSG::ERROR<<"PDFCutArray: Unable to open data file: "<<file<<endreq;
    clear();
    return;
  }

  std::string line;
  getline(dataFile, line);
  m_Narray= 0;
  m_Npdf= 0;
  int nArr= 0;
  for(; !dataFile.eof() && (!m_Narray || nArr<m_Narray || !m_Cuts);
        getline(dataFile, line) )
  {
    if( checkLine(line) ) continue;
    // used in CalLikelihoodManagerTool
    else if( checkField(line, "CUT TYPE") ) {
      if( checkField(line, "High Energy") ) 
        m_CutsFcn= new PDFHighEnergyCuts(tool, log);
      else m_CutsFcn= new PDFLowEnergyCuts(tool, log);
      continue;
    }
    else if( checkField(line, "#ARRAY") )
    {
      sscanf(line.data(), "#ARRAYS: %d #PDFs: %d %*s", &m_Narray, &m_Npdf);
      if( m_Narray>0 && m_Npdf>0 ) continue;
    }
    else if( checkField(line, "AXES:") )
    {
      m_Grid= new PDFGrid(dataFile, log);
      if( !(m_Grid[0]) ) 
        log<<MSG::ERROR<<"PDFCutArray: PDFGrid could not be created"<<endreq;
      else {
        m_LikelihoodFcn->setNinterpolationParameters(m_Grid->getNaxes());
        m_CutsFcn->setNinterpolationParameters(m_Grid->getNaxes());
        continue;
      }
    }
    else if( !m_Grid || !m_Narray || !m_CutsFcn )
      log<<MSG::ERROR<<"PDFCutArray: Missing PFDGrid, #Arrays,"
                       "or PDFFunction Cuts"<<endreq;

    else if( checkField(line, "CUTS:") )
    {
      m_Cuts= new PDFVertexArray(dataFile, this, m_CutsFcn, log);
      if( !(m_Cuts[0]) )
        log<<MSG::ERROR<<"PDFCutArray: Cuts PDFVertexArray was not created"
                       <<endreq;
      else continue;
    }
    else if( checkField(line, "PDF:") )
    {
      push_back(new PDFVertexArray(dataFile, this, m_LikelihoodFcn, log));
      if( !(back()[0]) )
        log<<MSG::ERROR<<"PDFCutArray: PDFVertexArray was not created"
                 <<endreq;
      else { ++nArr; continue; }
    }

    log<<MSG::ERROR<<"PDFCutArray: Incorrect data in file: "<<file<<endreq
             <<"                   At Array: "<<nArr<<endreq
             <<"                   Stopping on line: \""<<line<<"\""
             <<endreq;
    dataFile.close();
    clear();
    return;
  }
  dataFile.close();

  if( !m_Cuts ){
    log<<MSG::ERROR<<"PDFCutArray: No Cuts Defined"<<endreq;
    clear();
    return;
  }
  if( m_Narray!=nArr ){
    log<<MSG::ERROR<<"PDFCutArray: Missing PDFVertexArray"<<endreq;
    clear();
    return;
  }
}

void PDFCutArray::clear(void){
  if( m_Grid) { delete m_Grid; m_Grid= 0; }
  if( m_CutsFcn) { delete m_CutsFcn; m_CutsFcn= 0; }
  if( m_LikelihoodFcn) { delete m_LikelihoodFcn; m_LikelihoodFcn= 0; }
  if( m_Cuts) { delete m_Cuts; m_Cuts= 0; }
  PDFAnswer::clear();
}

bool PDFCutArray::evalEnergy(double eMin, double eMax,
                             const Event::CalCluster *cluster,
                             const Event::TkrVertex *vertex){
  setAnswer(0);
  m_LikelihoodFcn->setEvt(cluster, vertex);
  m_CutsFcn->setEvt(cluster, vertex);
  if( eMax<m_Grid->minTrialEnergy() || eMin>m_Grid->maxTrialEnergy() )
    return 0;
  if( eMin<m_Grid->minTrialEnergy() ) eMin= m_Grid->minTrialEnergy();
  if( eMax<m_Grid->getBinCenter(0, 4) ) eMax= m_Grid->getBinCenter(0, 4);
  else if( eMax>m_Grid->maxTrialEnergy() ) eMax= m_Grid->maxTrialEnergy();

  const double *fcn= m_LikelihoodFcn->getInterpolationParameters();
  if( m_Grid->initialise(fcn) ) return 0;
  return PDFAnswer::evalEnergy(eMin, eMax, cluster, vertex);
}

PDFAnswer::Status_t PDFCutArray::evalStatus(int i, PDFVertexArray *pdf){
  if( m_LikelihoodFcn->unPhysical(pdf->getEnergy()) ) return kUnphysical;
  double ans;
  m_CutsFcn->trialEnergy()= pdf->getEnergy();
  if( m_CutsFcn->eval((PDFParameters*) m_Cuts->at(i), &ans, true) )
    return kNoCuts;
  if( pdf!=((PDFVertexArray*) (at(int(ans)))) ) return kOutsideCut;
  return kInsideCut;
}

/******************************************************************************/
/***********************  PDFVertexArray **************************************/
/******************************************************************************/
PDFVertexArray::PDFVertexArray(std::ifstream &dataFile, PDFCutArray *tool,
                                            PDFFunction *fcn, MsgStream &log)
    :PDFAnswer(), m_Fcn(fcn), m_Tool(tool) {
  int pdf= 0;
  // Read in the parameters from the data file
  std::string line;
  getline(dataFile, line);
  int arrayLength= tool->getNpdf()/tool->getNarray();
  int nPar= m_Fcn->getNfunctionParameters();
  for( ; !dataFile.eof() && pdf<arrayLength; getline(dataFile, line) )
  {
    if( checkLine(line) ) continue;
    else if( checkField(line, "ARRAY") )
      log<<MSG::ERROR<<"PDFVertexArray: Reached next array without "
                       "filling this one"<<endreq;
    else if( checkField(line, "TABLE") )
    {
      push_back(new PDFParameters(dataFile, tool->getGrid(), nPar, log));
      ++pdf;
      if( !(back()[0]) )
        log<<MSG::ERROR<<"PDFVertexArray: PDF #"<<pdf-1
                 <<" did not build correctly"<<endreq;
      else continue;
    }

    log<<MSG::ERROR<<"PDFVertexArray: Incorrect data"<<endreq
                   <<"          At PDF: "<<pdf<<endreq
                   <<"          Stopping on line :\""<<line
                   <<"\""<<endreq;
    clear();
    return;
  }
}
  
bool PDFVertexArray::evalEnergy(double eMin, double eMax,
                          const Event::CalCluster*,
                          const Event::TkrVertex *vertex){
  setAnswer(0);
  m_MPV[0]= m_MPV[1]= 0.;
  m_status= PDFAnswer::kNotCalculated;
  int vZ0= getTkrPlane(vertex);

  if( m_Fcn->getMPV((PDFParameters *) at(vZ0), eMin, eMax, m_MPV) )
    return true;
  setAnswer(this);
  m_status= m_Tool?m_Tool->evalStatus(vZ0, this):kInsideCut;
  return false;
}

bool PDFVertexArray::evalError(Event::TkrVertex *vertex){
  m_FWHM[0]= m_FWHM[1]= 0.;
  return m_Fcn->getFWHM((PDFParameters*) at(getTkrPlane(vertex)),
                        m_MPV, m_FWHM);
}

/******************************************************************************/
/***********************  PDFAnswer *******************************************/
/******************************************************************************/
void PDFAnswer::clear(void){
  for(PDFVectItr pdf= begin(); pdf!=end(); ++pdf)
    if( *pdf ) delete *pdf;
  PDFVect::clear();
}

bool PDFAnswer::operator<(const PDFAnswer &a) const{ 
  if( getStatus()==a.getStatus() ) return getProbability()<a.getProbability();
  bool ans= getStatus()<a.getStatus();
  return ans;
}

bool PDFAnswer::evalEnergy(double eMin, double eMax,
                           const Event::CalCluster *cluster,
                           const Event::TkrVertex *vertex){
  m_Answer= 0;
  for(PDFVectItr pdf= begin(); pdf!=end(); ++pdf){
    PDFAnswer *ans= (PDFAnswer*)(*pdf);
    if(ans->evalEnergy(eMin, eMax, cluster, vertex)) continue;
    if( operator<(ans[0]) ) m_Answer= ans;
  }
  return !m_Answer;
}

bool PDFAnswer::checkLine(std::string &line) {
  return (line=="\n") || strspn(line.data(), " ")==strlen(line.data())
                      || checkField(line, "COMMENT");
} 
bool PDFAnswer::checkField(std::string &line, const char field[])
{ return line.size()>=strlen(field) && line.find(field)!=std::string::npos; }

int PDFAnswer::getAnswerID(void) const{
  if( !m_Answer || m_Answer==this ) return 0;
  
  int id= 1;
  for(PDFVectConstItr tool= begin();
      tool!=end() && ((PDFAnswer*) (*tool))!=m_Answer; ++tool, ++id);
  return id;
}

