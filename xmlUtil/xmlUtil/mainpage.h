// Mainpage for doxygen

/*! \mainpage package xml
<p>
This package provides services for manipulating certain generic xml
elements.  At this time these include elements for describing 
constants and arithmetic  (&lt;prim>, &lt;mul>, etc.) and those having 
to do with recording the genesis of an xml file (&lt;source>).</p>
<p>
This package also provides classes for representing and manipulating
identifier dictionaries and identifier conversions. An <b>identifier</b>
logically consists of a sequence of (name,value) pairs in which the
value is a non-negative integer.  An <b>identifier dictionary</b> is
a means of specifying constraints on a collection of Identifiers
in such a way that 
<ul>
  <li>The collection of identifiers form a tree</li>
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
<p>In addition to a test program, this facility includes stand-alone
programs to automatically generate special-purpose xml files from
an initial file describing an instrument.
</p>
<hr>
  \section requirements requirements
  \include requirements

*/
