//  $Header$

#include "dom/DOM_Element.hpp"
#include "dom/DOM_Document.hpp"
#include "dom/DOMString.hpp"

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/docMan/DocMan.h"
#include "xmlUtil/docMan/DocClient.h"

namespace xmlUtil {

  DocMan* DocMan::m_self = 0;

  DocMan::DocMan() {
    m_self = this;
    m_meFirst = 0;
    m_parser = new xml::XmlParser();
  }

  DocMan* DocMan::getPointer() {
    if (m_self == 0) {
      m_self = new DocMan();
    }
    return m_self;
  }

  bool DocMan::parse(const std::string& filename) {
    DOM_Document doc = m_parser->parse(filename.c_str());
    if (doc == DOM_Document()) return false;

    if (m_meFirst != 0) m_meFirst->handleChild(doc);

    DOM_Element root = doc.getDocumentElement();

    DOM_Element child = xml::Dom::getFirstChildElement(root);

    while (child != DOM_Element()) {
      DOMString eltName = child.getTagName();
      if (m_meFirst != 0) m_meFirst->handleChild(child);

      ClientList* list = findList(std::string(xml::Dom::transToChar(eltName)));
      if (list) list->invoke(child);

      child = xml::Dom::getSiblingElement(child);
    }

    // Clean up call for privileged client
    if (m_meFirst != 0) m_meFirst->handleChild(DOM_Node());

    m_parser->reset();
    return true;
  }

  bool DocMan::regClient(const std::string& eltName, DocClient* client) {
    ClientList* curList = findList(eltName);
    if (curList == 0) {
      curList = new ClientList(eltName);
      m_lists.push_back(curList);
    }
    return curList->add(client);
  }

  bool DocMan::regMeFirst(DocClient* client) {
    if (m_meFirst != 0)       return false;

    m_meFirst = client;
    return true;
  }

  DocMan::ClientList* DocMan::findList(const std::string& eltName) {
    unsigned int ix;

    for (ix = 0; ix < m_lists.size(); ix++) {
      ClientList* cur = m_lists[ix];
      if (eltName.compare(cur->getName()) == 0) return cur;
    }
    return 0;
  }

  bool DocMan::ClientList::add(DocClient* client, bool front) {
    ClientsIt it = m_clients.begin();
    const std::string& clientName = client->getName();
    while (it != m_clients.end() ) {
      // can't register for the same element type more than once
      if (clientName.compare((*it)->getName()) == 0 ) return 0;
      ++it;
    }
    if (front) m_clients.insert(m_clients.begin(), client);
    else m_clients.push_back(client);

    return true;
  }

  // This is just here for completeness.  Don't need a real 
  // implementation yet because no one calls it
  bool DocMan::ClientList::remove(DocClient* client) {
    return true;
  }

  void DocMan::ClientList::invoke(DOM_Element elt) {
    ClientsIt it = m_clients.begin();
    while (it != m_clients.end() ) {
      (*it)->handleChild(elt);
      ++it;
    }
  }

  // Still need destructors
  DocMan::ClientList::~ClientList() {
    m_clients.clear();
  }

  DocMan::~DocMan() {
    ListsIt it = m_lists.begin();
    while (it != m_lists.end()) {
      delete (*it);
      ++it;
    }
    delete m_parser;
  }
}
