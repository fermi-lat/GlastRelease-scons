
#ifndef IMEDIA_H
#define IMEDIA_H

/** @class IMedia
* @brief Abstract class that defines callbacks from a traversal of the
* detector materials
*
* $Header$
*/
class IMedia {
public:

    typedef std::vector<double> DoubleVector;
    typedef std::vector<unsigned int>UintVector;
    typedef std::vector<std::string>StringVector;

    enum MediaType{ Element, Material};
    enum CompositeType{ Natoms, Fract};

    /** Create a new simple material
        @param name of the material
        @param type type of the material
	@param params Parameters vector 
	                    - (z,a) for element
			    - (z,a,density) for material
	@param symbol The symbol (if applicable)
        */
    virtual void addMaterial(std::string name, 
			     MediaType type,
			     const DoubleVector& params, std::string symbol)=0;

    /** Create a new composite material out of a list of materials name
        @param name of the material
        @param type type of the composite (natoms or fractional)
	@param density The density of the composite
	@param components The list of components names
	@param qty The vector of quantities (it is interpreted as int or 
	           double depending on type)
        */
    virtual void addComposite(std::string name, 
			      CompositeType type,
			      double density,
			      const StringVector& components, 
			      const DoubleVector& qty)=0;
    

    virtual ~IMedia(){}

};
#endif
