// Mainpage for doxygen

/*! \mainpage package xmlUtil
<p>
This package provides services for manipulating certain generic xml
elements.  At this time these include elements for describing 
constants and arithmetic  (<prim>, <mul>, etc.) and those having 
to do with recording the genesis of an xml file (<source>).</p>
<p>The <b>DocMan</b> class and the derived class <b>GDDDocMan</b>
facilitate sharing of a single XML document among several clients.
Clients can sign up to be called back when an element of a particular
type is encountered.  DocMan itself is responsible for parsing the
file and managing the resources this requires. GDDDocMan "knows"
about constants; it provides the additional service of evaluating
and substituting for constants before any client handlers are invoked.</p>
<p>
This package also provides classes for representing and manipulating
identifier dictionaries and identifier conversions. An <b>identifier</b>
logically consists of a sequence of (name,value) pairs in which the
value is a non-negative integer.  An <b>identifier dictionary</b> is
a means of specifying constraints on a collection of Identifiers
in such a way that 
<ul>
  <li>The collection of identifiers forms a tree</li>
  <li>Given the value sequence from an identifier, it is possible
      to reconstruct the (unique) corresponding name sequence.</li>
</ul></p>
<p>An <b>identifier conversion</b> is a list of operations defined on
a subcollection of identifiers from a dictionary, transforming them 
to a new collection of identifiers (not necessarily belonging to the same 
dictionary).  An <b>id converter</b> is a group of compatible identifier
conversions. </p>

<p>
These elements (for arithmetic, identifiers, etc.) are described  
in gdd.dtd, the dtd describing generic
geometry, hence that file and the xml files using it are also kept
in this package.  When xml namespaces become usable with a c++
parser it would be preferable to split off the utility elements
into their own namespace(s). 
</p>
<p>In addition to several test program, this facility includes stand-alone
programs to automatically generate special-purpose xml files from
an initial file describing an instrument.
</p>
<p>See release.notes for recent activity.</p>
<hr>
  \section requirements requirements
  \include requirements

*/
