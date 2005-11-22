/**@file XprsnTree.h
@brief Contains class definitions for implementing a very basic Decision Tree
@author T. Usher
$Header$
*/

#ifndef XprsnTree_h
#define XprsnTree_h

#include <iostream>

/** @class Exception 
    @brief hold a string
*/
class XTENexception : public std::exception
{
public: 
    XTENexception(std::string error):m_what(error){}
    ~XTENexception() throw() {;}
    virtual const char *what( ) const  throw() { return m_what.c_str();} 
    std::string m_what;
};

/** @class IXTEprsnNode
@brief Interface definition for an XTExprsnNode 
*/
class IXTExprsnNode 
{
public:
    // Return the value contained here (to be overidden at higher levels)
    virtual const void* operator()() const = 0;

    // All derived objects will be "named" so return their name here
    virtual const std::string& getName() const = 0;
    
    // Basic print full expresssions to stream method
    virtual void print(std::ostream& out=std::cout, bool first=true) const = 0;

    // This is for printing name and value
    virtual void printExp(std::ostream& out=std::cout, bool first=true) const = 0;
};


/** @class XTExprsnNode
@brief A generic implementation of an "Expression Node". This class defines the 
       actual type of the node and implements arithmetic or logical functions
       which can be performed there. 
       The interface class requires the value to be returned as a void*... which
       needs to be "reinterpret_cast" back to type T by the caller. 

       The Template parameter deterimes the type of output we are. 
       In actual use this will be either bool or double
*/
template <class T> class XTExprsnNode : virtual public IXTExprsnNode
{
public:

    // The basic constructor 
    XTExprsnNode(const std::string& name,
                 IXTExprsnNode&     left, 
                 IXTExprsnNode&     right) 
                  : m_name(name), m_value(new T), m_left(left), m_right(right) 
    {
        // Set the function pointer to the type of operation we define
        if      (name == "+")  mathOp = &XTExprsnNode<T>::operator+;  // Arithmetic addition
        else if (name == "-")  mathOp = &XTExprsnNode<T>::operator-;  // Arithmetic subraction
        else if (name == "*")  mathOp = &XTExprsnNode<T>::operator*;  // Arithmetic multiplication
        else if (name == "/")  mathOp = &XTExprsnNode<T>::operator/;  // Arithmetic division
        else if (name == "|")  mathOp = &XTExprsnNode<T>::operator||; // Logical "OR" (support IM mode)
        else if (name == "||") mathOp = &XTExprsnNode<T>::operator||; // Logical "OR" 
        else if (name == "&")  mathOp = &XTExprsnNode<T>::operator&&; // Logical "AND" (support IM mode)
        else if (name == "&&") mathOp = &XTExprsnNode<T>::operator&&; // Logical "AND"
        else if (name == ">")  mathOp = &XTExprsnNode<T>::operator>;  // Logical compare >
        else if (name == ">=") mathOp = &XTExprsnNode<T>::operator>=; // Logical compare >=
        else if (name == "<")  mathOp = &XTExprsnNode<T>::operator<;  // Logical compare <
        else if (name == "<=") mathOp = &XTExprsnNode<T>::operator<=; // Logical compare <=
        else if (name == "==") mathOp = &XTExprsnNode<T>::operator==; // Logical compare ==
        else if (name == "!=") mathOp = &XTExprsnNode<T>::operator!=; // Logical compare != 
        else if (name == "^")  mathOp = &XTExprsnNode<T>::pow;        // Raise x to power y 
        else throw XTENexception("XTExprsnNode: Invalid math operation requested");
    }
    virtual ~XTExprsnNode() {delete m_value;}

    // Overide the generic "retrieve" operator from the base class
    virtual const void* operator()() const     {return (this->*mathOp)(m_right);}
    virtual const std::string& getName() const {return m_name;}

    // Provide mechanism for outputting to a standard stream
    virtual void print(std::ostream& out=std::cout, bool first=true) const
    {
        // Go down the left branch first
        m_left.print(out,false);

        // Print out name
        out << " " << m_name.data() << " ";

        // Now the right branch
        m_right.print(out,false);

        if (first) out << std::endl;
    }
    virtual void printExp(std::ostream& out=std::cout, bool first=true) const
    {
        out << "(";
        m_left.printExp(out,false);
        out << ")" << m_name.data() << "(";
        m_right.printExp(out,false);
        out << ")=" << (*m_value);
        if (first) out << std::endl;
    }

private:

    // The allowed operations:
    const T* operator+ (IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left +  rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) +  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator- (IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left - rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) -  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator* (IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left * rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) *  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator/ (IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left / rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) /  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator||(IXTExprsnNode& rhs) const 
    {
        *m_value = *(reinterpret_cast<const bool*>(m_left())) ||  *(reinterpret_cast<const bool*>(rhs()));
        return m_value;
    }
    const T* operator&&(IXTExprsnNode& rhs) const 
    {
        *m_value = *(reinterpret_cast<const bool*>(m_left())) &&  *(reinterpret_cast<const bool*>(rhs()));
        return m_value;
    }
    const T* operator> (IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left > rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) >  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator>=(IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left >= rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) >=  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator< (IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left < rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) <  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator<=(IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left <= rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) <=  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator==(IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left == rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) ==  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* operator!=(IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = left != rght;
        //*m_value = *(reinterpret_cast<const double*>(m_left())) !=  *(reinterpret_cast<const double*>(rhs()));
        return m_value;
    }
    const T* pow(IXTExprsnNode& rhs) const 
    {
        double left = *(reinterpret_cast<const double*>(m_left()));
        double rght = *(reinterpret_cast<const double*>(rhs()));
        *m_value = std::pow(left,rght);
        //*m_value = std::pow(*(reinterpret_cast<const double*>(m_left())),*(reinterpret_cast<const double*>(rhs())));
        return m_value;
    }
    
    // Function pointer to the particular operation defined for this node
    const T* (XTExprsnNode<T>::*mathOp)(IXTExprsnNode& rhs) const;

    // Our name and value
    std::string m_name;     // The name assigned to this node
    T*          m_value;    // The "value" of this node 

    // The left and right branches beneath the node 
    IXTExprsnNode& m_left;  // Left branch to follow
    IXTExprsnNode& m_right; // Right branch to follow

    // This is for making a fancy output... I'm not necessarily proud of it...
    //template <class T> friend std::ostream& operator <<(std::ostream& stream, const XTExprsnNode<T>& node);
    friend std::ostream& operator <<(std::ostream& stream, const XTExprsnNode<T>& node);
};

// Override the ostream << operator. Keeping here as an example....
template <class T> std::ostream& operator<<(std::ostream& stream, const XTExprsnNode<T>& node)
{
    node.printExp(out,false);
    return stream;
}

/** @class XTExprsnNode
@brief A generic implementation of an "Expression Node". This class defines the 
       actual type of the node and implements arithmetic or logical functions
       which can be performed there. 
       The interface class requires the value to be returned as a void*... which
       needs to be "reinterpret_cast" back to type T by the caller. 

       The Template parameter deterimes the type of output we are. 
       In actual use this will be either bool or double
*/
template <class T> class XTIfElseNode : virtual public IXTExprsnNode
{
public:

    // The basic constructor 
    XTIfElseNode(const std::string& name,
                 IXTExprsnNode&     condNode,
                 IXTExprsnNode&     ifNode, 
                 IXTExprsnNode&     elseNode) 
                  : m_name(name), 
                    m_value(new T), 
                    m_condNode(condNode),
                    m_ifNode(ifNode), 
                    m_elseNode(elseNode) 
    { }
    virtual ~XTIfElseNode() {delete m_value;}

    // Overide the generic "retrieve" operator from the base class
    virtual const void* operator()() const 
    {
        bool condition = *(reinterpret_cast<const bool*>(m_condNode()));
        if (condition) *m_value = *(reinterpret_cast<const T*>(m_ifNode()));
        else           *m_value = *(reinterpret_cast<const T*>(m_elseNode()));

        return m_value;
    }
    virtual const std::string& getName() const {return m_name;}

    // Provide mechanism for outputting to a standard stream
    virtual void print(std::ostream& out=std::cout, bool first=true) const
    {
        out << " ifelse(";
        m_condNode.print(out,false);
        out << ",";
        m_ifNode.print(out,false);
        out << ",";
        m_elseNode.print(out,false);
        out << ") ";

        if (first) out << std::endl;
    }

    virtual void printExp(std::ostream& out=std::cout, bool first=true) const
    {
        out << "Node: " << m_name << " ";
        m_condNode.printExp(out, first);
        if (first) out << std::endl;

        bool condition = *(reinterpret_cast<const bool*>(m_condNode()));
        if (condition) m_ifNode.printExp(out,true);
        else           m_elseNode.printExp(out,true);
    }

private:

    // Our name and value
    std::string m_name;     // The name assigned to this node
    T*          m_value;    // The "value" of this node 

    // The left and right branches beneath the node 
    IXTExprsnNode& m_condNode;  // Condition to test for
    IXTExprsnNode& m_ifNode;    // what to return if condition is true
    IXTExprsnNode& m_elseNode;  // what to return if condition is false

    // This is for making a fancy output... I'm not necessarily proud of it...
    //template <class T> friend std::ostream& operator <<(std::ostream& stream, const IXTExprsnNode& node);
    //template <class T> friend std::ostream& operator <<(std::ostream& stream, const XTIfElseNode<T>& node);
};

//template <class T> std::ostream& operator <<(std::ostream& stream, const XTExprsnNode<T>& node)
//template <class T> std::ostream& operator<<(std::ostream& stream, const XTIfElseNode<T>& node)
//{
//    const XTExprsnNode<bool>& condNode = dynamic_cast<const XTExprsnNode<bool>&>(node.m_condNode);
//    stream << "Node: " << node.m_name << " ";
//    stream << condNode << std::endl;
//    if (*node.m_value)
//    {
//        const XTExprsnNode<T>& ifNode = dynamic_cast<const XTExprsnNode<T>&>(node.m_ifNode);
//        stream << ifNode << std::endl;
//    }
//    else
//    {
//        const XTExprsnNode<T>& elseNode = dynamic_cast<const XTExprsnNode<T>&>(node.m_elseNode);
//        stream << elseNode << std::endl;
//    }
//    return stream;
//}


/** @class XTEprsnValue
@brief A base class for defining the nodes in a Decision Tree
       Defines basic data members (the name and its "value"), 
       a method for retrieving the value and, finally, a method
       to print it. 
*/
template <class T> class XTExprsnValue : virtual public IXTExprsnNode
{
public:
    XTExprsnValue<T>(const std::string& name, T* value) : m_name(name), m_value(value) {}
    virtual ~XTExprsnValue() {}

    // "My" functions
    const T getValue() const {return *m_value;}

    // Return the value contained here (to be overidden at higher levels)
    virtual const void* operator()() const {return m_value;}
    virtual const std::string& getName() const {return m_name;}
    
    // Basic print to stream method
    virtual void print(std::ostream& out=std::cout, bool first=true) const
    {
        out << m_name.data();
        if (first) out << std::endl;
    }

    virtual void printExp(std::ostream& out=std::cout, bool first=true) const
    {
        out << m_name.data() << "=" << *m_value;
        if (first) out << std::endl;
    }

protected:
    std::string m_name;     // The name assigned to this node
    T*          m_value;    // The "value" of this node 

    // This is for making a fancy output... I'm not necessarily proud of it...
    //template <class T> friend std::ostream& operator <<(std::ostream& stream, const IXTExprsnNode& node);
    friend std::ostream& operator <<(std::ostream& stream, const IXTExprsnNode& node);
};

/** @class XTEprsnTupleVal
@brief A base class for defining the nodes in a Decision Tree
       Defines basic data members (the name and its "value"), 
       a method for retrieving the value and, finally, a method
       to print it. 
*/
template <class T> class XTExprsnTupleVal : virtual public IXTExprsnNode
{
public:
    XTExprsnTupleVal<T>(const std::string& name, T* value) : m_name(name), m_value(value) {}
    virtual ~XTExprsnTupleVal() {}

    // Return the value contained here (to be overidden at higher levels)
    virtual const void* operator()() const {return (*m_value)();}
    virtual const std::string& getName() const {return m_name;}
    
    // Basic print to stream method
    virtual void print(std::ostream& out=std::cout, bool first=true) const
    {
        out << m_name.data();
        if (first) out << std::endl;
    }

    virtual void printExp(std::ostream& out=std::cout, bool first=true) const
    {
        out << m_name.data();
        out << "=" << *(*m_value)();
        if (first) out << std::endl;
    }

protected:
    std::string m_name;     // The name assigned to this node
    T*          m_value;    // The "value" of this node 
};


/** @class XTExprsnFunction
@brief 
*/
#include <cmath>

template <class T> class XTExprsnFunction : virtual public IXTExprsnNode
{
public:

    // The basic constructor 
    XTExprsnFunction(const std::string& name,
                     IXTExprsnNode&     left, 
                     IXTExprsnNode&     right) 
                    : m_name(name), m_value(new T), m_left(left), m_right(right) 
    {
        // Set the function pointer to the type of operation we define
        if      (name == "^")     mathOp = &XTExprsnFunction<T>::pow;   // x to power y
        else if (name == "log")   mathOp = &XTExprsnFunction<T>::log10; // log base 10
        else if (name == "loge")  mathOp = &XTExprsnFunction<T>::loge;  // natural log
        else if (name == "sqrt")  mathOp = &XTExprsnFunction<T>::sqrt;  // square root
        else if (name == "sin")   mathOp = &XTExprsnFunction<T>::sin;   // sine
        else if (name == "cos")   mathOp = &XTExprsnFunction<T>::cos;   // cosine
        else if (name == "tan")   mathOp = &XTExprsnFunction<T>::tan;   // tangent
        else if (name == "exp")   mathOp = &XTExprsnFunction<T>::exp;   // e raised to power x
        else if (name == "abs")   mathOp = &XTExprsnFunction<T>::abs;   // absolute value
        else if (name == "min")   mathOp = &XTExprsnFunction<T>::min;   // smallest
        else if (name == "max")   mathOp = &XTExprsnFunction<T>::max;   // largest
        else throw XTENexception("XTExprsnFunction: Invalid function requested");
    }
    virtual ~XTExprsnFunction() {delete m_value;}

    // Overide the generic "retrieve" operator from the base class
    virtual const void* operator()() const     {return (this->*mathOp)(m_left);}
    virtual const std::string& getName() const {return m_name;}

    // Provide mechanism for outputting to a standard stream
    virtual void print(std::ostream& out=std::cout, bool first=true) const
    {
        // Print out name
        out << " " << m_name.data() << "(";

        // Go down the left branch first
        m_left.print(out,false);

        // The right operand?
        if (m_right.getName() != "")
        {
            out << ",";
            // Now the right branch
            m_right.print(out,false);
        }

        out <<") ";

        if (first) out << std::endl;
    }
    virtual void printExp(std::ostream& out=std::cout, bool first=true) const
    {
        out << m_name.data() << "()";
        //out << "=" << (*m_value);
        if (first) out << std::endl;
    }

private:

    // The allowed operations:
    const T* pow(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        double rght = *(reinterpret_cast<const double*>(m_right()));
        *m_value = std::pow(left,rght);
        //*m_value = std::pow(*(reinterpret_cast<const double*>(arg())),*(reinterpret_cast<const double*>(m_right())));
        return m_value;
    }
    const T* log10(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        *m_value = std::log(left);
        //*m_value = std::log10(*(reinterpret_cast<const double*>(arg())));
        //*m_value = std::log(*(reinterpret_cast<const double*>(arg())));
        return m_value;
    }
    const T* loge(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        *m_value = std::log(left);
        //*m_value = std::log(*(reinterpret_cast<const double*>(arg())));
        return m_value;
    }
    const T* sqrt(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        *m_value = std::sqrt(left);
        //*m_value = std::sqrt(*(reinterpret_cast<const double*>(arg())));
        return m_value;
    }
    const T* sin(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        *m_value = std::sin(left);
        //*m_value = std::sin(*(reinterpret_cast<const double*>(arg())));
        return m_value;
    }
    const T* cos(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        *m_value = std::cos(left);
        //*m_value = std::cos(*(reinterpret_cast<const double*>(arg())));
        return m_value;
    }
    const T* tan(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        *m_value = std::tan(left);
        //*m_value = std::tan(*(reinterpret_cast<const double*>(arg())));
        return m_value;
    }
    const T* exp(IXTExprsnNode& arg) const 
    {
        double left = *(reinterpret_cast<const double*>(arg()));
        *m_value = std::exp(left);
        //*m_value = std::exp(*(reinterpret_cast<const double*>(arg())));
        return m_value;
    }
    const T* abs(IXTExprsnNode& arg) const 
    {
        double temp = *(reinterpret_cast<const double*>(arg()));
        *m_value = temp > 0 ? temp : -temp;
        return m_value;
    }
    const T* min(IXTExprsnNode& arg) const 
    {
        double left  = *(reinterpret_cast<const double*>(arg()));
        double right = *(reinterpret_cast<const double*>(m_right()));
        *m_value = left < right ? left : right;
        return m_value;
    }
    const T* max(IXTExprsnNode& arg) const 
    {
        double left  = *(reinterpret_cast<const double*>(arg()));
        double right = *(reinterpret_cast<const double*>(m_right()));
        *m_value = left > right ? left : right;
        return m_value;
    }
    
    // Function pointer to the particular operation defined for this node
    const T* (XTExprsnFunction<T>::*mathOp)(IXTExprsnNode& rhs) const;

    // Our name and value
    std::string m_name;     // The name assigned to this node
    T*          m_value;    // The "value" of this node 

    // The left and right branches beneath the node 
    IXTExprsnNode& m_left;  // Left branch to follow
    IXTExprsnNode& m_right; // Right branch to follow
};

#endif
