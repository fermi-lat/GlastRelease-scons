// $Header$

#ifndef XMLUTIL_IDOPFACTORY_H
#define XMLUTIL_IDOPFACTORY_H


namespace xmlUtil {

  class IdOperation;

  //! Factory which knows how to make any kind of IdOperation.
  /*! It's a singleton because there is no need for more than one */
  class IdOpFactory {
  public:

    //! Return pointer to the one and only instance, creating it first
    //! if necessary
    static IdOpFactory* getIdOpFactory() {
      if (!me) me = new IdOpFactory();
      return me;
    }

    //! Depending on which DOM_Element we're passed, instantiate
    //! an object from the appropriate IdOperation derived class.
    IdOperation *makeIdOp(DOM_Element);

  private:
    //! Private constructor for singleton
    IdOpFactory(){};

    static IdOpFactory* me = 0;
  };  // End IdOpFactory class
}   // end xmlUtil namespace
