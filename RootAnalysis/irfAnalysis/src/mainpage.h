//$Header$
// (Special "header" just for doxygen)

/*! @mainpage  package irfAnalysis

@authors T. Burnett         <TBurnett@u.washington.edu>, 
         Jim Chiang         <JChiang@slac.stanford.edu>, 
         T. Hansl-Kozanecka <hansl@slac.stanford.edu>

Builds several small applications to display instrument response
 - aeffTable
 - aeff
 - energy_fit
 - psf_fit

The paths to the input and output files are defined in the requirements. 
Output data (psf_friend.root, aeffTable.root, etc) are saved in subdirectory 
<tt>data</tt>. 
  
<hr>
<h3> Classes common to applications  </h3>
<ul>
  <li>IRF Base class <ul> 
      <li>  opens the input file, 
      <li>  defines energy and angle grid and the corresponding histogram names,
      <li>  and provides utility for canvas layout. </ul>
  <li>PSF Interface for applications <ul>
      <li> implements opening of the input file and 
           creates the friend ROOT file with derived quantities;
      <li> several methods to fill and draw plots.            
    </ul>
</ul>
<hr>
<h3> aeffTable  </h3>
AeffTable 
<hr>
<h3> aeff       </h3>
<h3> energy_fit </h3>
<h3> psf_fit    </h3>

<hr>
\section requirements requirements
\include requirements
<hr>
    @todo Several numbers, which are used for the generation of events and 
          are needed to calculate the effective area, 
          are hard coded, but should be read from a bookkeeping file: 
          energy range, number of events generated, target area used. 
          

*/
