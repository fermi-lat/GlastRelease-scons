// $Header$

#ifndef XMLUTIL_IDCONVERTER_H
#define XMLUTIL_IDCONVERTER_H


namespace xmlUtil {
  class IdDict;
  class Identifier;
  class IdConversion;

  class IdConverter {
  public:
    IdConverter(DOM_Element elt);
    ~IdConverter();

    //! Return image of the input dictionary under the conversions
    //! defined in the IdConverter object.
    /*! \arg \b inDict: the input dictionary
        \arg \b status: report on success of conversion, validity of
             output dictionary (and hence of IdConverter object).
        Note the collection of conversions making up a converter must
        act on disjoint domains, which is identical to the restriction
        that no conversion have a path which is a subset of another
        converter's path.
    */
    IdDict *convert(const IdDict const *inDict, unsigned *status) const;

    //! Return the result of converting a single NamedId
    NamedId *convert(const NamedId const *in) const;

    //! Return the result of converting a single Identifier.
    /*! Can only do this if there is an input IdDict so that
        the IdDict can be asked for the NamedId corresponding
        to the Identifier */
    Identifier *convert(const Identifier const *inIdent) const;

    // Check that conversions comprising the converter have disjoint
    // domains.
    bool isConsistent();

  private:
    // Should inputDict be optional or required?
    IdDict *inputDict; /*< Input IdDict.   */

    // Might choose a different container for the conversions, for
    // example a tree, ordered by input path on which a conversion
    // operates
    typedef std:vector<IdConversion*>  Conversions; 

    Conversions conv;  /*< collection of conversion operations making up the 
                         IdConverter */

  };  // end class IdConverter

}  // end namespace




#endif
