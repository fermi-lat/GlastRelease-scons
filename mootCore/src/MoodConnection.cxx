// $Header$

#include <string>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "mootCore/MoodConnection.h"
#include "rdbModel/Db/Connection.h"
#include "rdbModel/Db/MysqlConnection.h"
#include "rdbModel/Db/ResultHandle.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Tables/Column.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/RdbException.h"
#include "rdbModel/Management/Manager.h"
#include "rdbModel/Management/XercesBuilder.h"
#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"
#include "facilities/Util.h"
#include "facilities/commonUtilities.h"
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>

namespace MOOT {

  MoodConnection::MoodConnection(bool wrt, bool dbg) : m_rdb(0), m_write(wrt),
                                                       m_dbg(dbg) {
    using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
    using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
    using xmlBase::Dom;
    using facilities::Util;

    // Always look up host, port, dbname from myHost.xml in archive,
    // for both read and write
    std::string hostfile("$(MOOT_ARCHIVE)/myHost.xml");
    xmlBase::XmlParser parser;
    XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument* doc;
    parser.doSchema(true);

    try {
      doc = parser.parse(hostfile.c_str());
    }
    catch (std::exception ex) {
      std::cerr << ex.what();
      std::cerr.flush();
      throw ex;
    }
    if (doc == NULL) {
      throw std::runtime_error("Archive does not have proper myHost.xml file");
    }
    DOMElement* connElt = doc->getDocumentElement();
    std::string dbstring = Dom::getAttribute(connElt, "dbname");
    m_nick = Dom::getAttribute(connElt, "nickname");
    DOMElement* masterElt = Dom::getFirstChildElement(connElt);
    std::string hostStr;
    std::string portStr;
    if (!wrt) {  // try slave if there is one
      DOMElement* slaveElt = Dom::getSiblingElement(masterElt);
      if (slaveElt) {
        hostStr = Dom::getAttribute(slaveElt, "host");
        portStr = Dom::getAttribute(slaveElt, "port");
        bool ok= open(hostStr.c_str(), 
                      Util::stringToInt(portStr), dbstring, wrt);
        if (ok) return;
        else {
          std::cerr << "Failed to open connection to MOOT db" << std::endl;
          std::cerr.flush();
          exit(1);
        }
      }
    }
    else {   // asked for write. Must be real archive, not a mirror
      if (Dom::getAttribute(connElt, "write") != "true") { 
        // tilt!
        throw std::runtime_error("Not allowed to write to slave dbs");
      }
    }
    // If we got this far, must use master host, port
    hostStr = Dom::getAttribute(masterElt, "host");
    portStr = Dom::getAttribute(masterElt, "port");
    if (!(open(hostStr.c_str(), Util::stringToInt(portStr), 
               dbstring, wrt)) )                               {
      std::cerr << "Failed to open connection to MOOT db" << std::endl;
      std::cerr.flush();
      exit(1);
    }
    return;
  }


  MoodConnection::~MoodConnection() {
    if (m_rdb) {
      close();
    }
  }

  void MoodConnection::close() {
    if (m_rdb) {
      rdbModel::Connection* conn = m_rdb->getConnection();
      conn->close();
      delete conn;
      delete m_rdb;
      m_rdb = 0;
    }
  }

  bool MoodConnection::open(const char* host, int port, 
                            std::string& dbname, bool wrt) {
    using facilities::commonUtilities;
    static const char* acct = "glastreader";
    // If wrt is true, expect to get user, password from .my.cnf
    // Else use glastreader

    if (m_rdb) {
      if (wrt == m_write) {
        std::cout << 
          "MoodConnection::oldOpen  Warning: already have open connection" 
                  << std::endl;
        return true;
      }
      else {    // close old connection first
        close();
      }
    }
    
    m_write = wrt;
    std::ostream* infoOut = 0;   // goes to std::cout
    rdbModel::MysqlConnection* conn = 0;
    if (!m_dbg) {
      infoOut = new std::ofstream("/dev/null");
    }
    conn = new rdbModel::MysqlConnection(infoOut);
    if (!conn) return false;
    bool ok = conn->init();
    if (!ok) {
      delete conn;
      return false;
    }
    std::string cnf("$(HOME)/.my.cnf");
    if (m_write) {
      facilities::Util::expandEnvVar(&cnf);
      ok = conn->setOption(rdbModel::DBreadDefaultFile, 
                             cnf.c_str());
      if (ok) ok = conn->setOption(rdbModel::DBreadDefaultGroup, "MOOT_WRITE");
      if (ok) ok =  conn->open(host, port, NULL, NULL, dbname.c_str());
    }
    else { // use glastreader
      ok = conn->open(host, port, acct, acct, dbname.c_str());
    }
    if (!ok) {
      delete conn;
      m_rdb = 0;
      return false;
    }
    m_rdb = new rdbModel::Rdb;
    m_rdb->setConnection(conn);
    // look for xml file describing db
    // If MOOT_XML has a value, use it.  Otherwise, all we need is our 
    // own xml dir.
    std::string dbSchema = commonUtilities::getEnvironment("MOOT_XML");
    if (dbSchema.size() == 0) {
      dbSchema = commonUtilities::getXmlPath("mootCore");
    }
#ifdef WIN32
    std::string slash("\\");  // need to escape properly
#else
    std::string slash("/");    
#endif
    dbSchema += slash + dbname + ".xml";

    // expandEnvVar throws Untranslatable exception if 
    try {
      facilities::Util::expandEnvVar(&dbSchema);
    }
    catch (facilities::Untranslatable ex) { 
      std::cerr << "Missing definition for environment variable MOOT_XML"
                << std::endl;
      std::cerr.flush();
      delete conn;
      m_rdb = 0;
      throw ex;
    }

    int status = compareSchema(dbSchema);
    if (status != 1) {
      delete conn;
      m_rdb = 0;
      throw 
        std::invalid_argument("Bad db description file or lost connection");
    }
    return true;
  }


  int MoodConnection::compareSchema(const std::string& schema) {
    using namespace rdbModel;

    // Maybe first check if file exists?  E.g., try opening for read

    // good errCode is 0
    rdbModel::Builder* b = new rdbModel::XercesBuilder;
    int errCode = m_rdb->build(schema, b);

    if (errCode) {
      std::cerr << "Error in database description file " << schema 
                << std::endl;
      std::cerr << "Parse failed with error " << errCode << std::endl;
      return -1; // RETBadCnfFile;
    }

    rdbModel::MATCH match = m_rdb->getConnection()->matchSchema(m_rdb, false);

    switch (match) {
    case rdbModel::MATCHequivalent:
    case rdbModel::MATCHcompatible:
      if (m_dbg) {
        std::cout << "XML schema and MySQL database are compatible" 
                  << std::endl;
      }
      return 1; // RETOk;
    case rdbModel::MATCHfail:
      std::cerr << "XML schema and MySQL database are NOT compatible" 
                << std::endl;
      return -1; // RETBadCnfFile;
    case rdbModel::MATCHnoConnection:
      std::cerr << "Connection failed while attempting match" << std::endl;

      return -2; // RETNoConnect;
    }
    return -3; // RETBadValue;
  }
}
