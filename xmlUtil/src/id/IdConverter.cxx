// $Header$

#include "xmlUtil/id/IdConversion.h"
#include "xmlUtil/id/IdConverter.h"
#include "xmlUtil/id/IdDictMan.h"
#include "xmlUtil/id/IdDict.h"
#include "dom/DOM_Element.hpp"
#include "xml/Dom.h"

namespace xmlUtil {

  IdDictMan * IdConverter::dictMan = 0;

  IdConverter::IdConverter(DOM_Element elt) : consistent(UNKNOWN),
  sorted(NO) {
    if (!dictMan) dictMan = IdDictMan::getPointer();

    // An idConverter element may optionally have as first child
    // a constants element, but we're assuming that the XML document
    // has already had constants evaluation and substitution done
    // on it, if necessary.
    // Assume that any constant evaluations and substitutions have
    // already been dealt with 
    DOM_Element child = xml::Dom::getFirstChildElement(elt);

    if ( (child.getTagName()).equals(DOMString("constants")) ) { // move on
      child = xml::Dom::getSiblingElement(child);
    }
    while ( (child.getTagName()).equals(DOMString("idDict")) ) {
      IdDict* dict = new IdDict(child);
      dictMan->registerDict(dict);
      child = xml::Dom::getSiblingElement(child);
    }

    // Now that all id dictionaries we need should have been defined,
    // handle attributes
    std::string dictName = xml::Dom::getAttribute(elt, "fromDict");
    if (!(dictName.compare("")) ) {
      inputDictName = new std::string(dictName);
      inputDict = dictMan->findDict(dictName);
    }
    dictName = xml::Dom::getAttribute(elt, "toDict");
    if (!(dictName.compare("")) ) {
      outputDictName = new std::string(dictName);
      outputDict = dictMan->findDict(dictName);
    }

    // All remaining elements should be idConv
    while (child != DOM_Element() ) {
      IdConversion *conv = new IdConversion(child);
      convCol.push_back(conv);
      child = xml::Dom::getSiblingElement(child);
    }

  }

  IdConverter::~IdConverter() {
    convCol.clear();
  }

  //   *** TO DO ****  Put in real code. Use std::sort algorithm  ??
  void IdConverter::sortConvs() {
    if (sorted != YES) {
      ;  // sort them
      sorted = YES;
    }
  }

  //   *** TO DO ****  Put in some real code!
  bool IdConverter::isConsistent() {
    if (consistent == UNKNOWN) {
      sortConvs();
      // Check whether a path is a subsequence of some other path
      consistent = YES;
    }
    return consistent;
  }

  //   *** TO DO ***    (for now just make a copy)
  NamedId *IdConverter::convert(const NamedId *in) const {
    return new NamedId(*in);
  }

  //   *** TO DO ***    (for now just check for inDict; make a copy)

  Identifier *IdConverter::convert(const Identifier *inIdent) const {
    if (!inputDict) return 0;
    return new Identifier(*inIdent);
  }
}

