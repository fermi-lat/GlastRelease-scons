#ifndef ROOT_GCRSELECT_H
#define ROOT_GCRSELECT_H

#include "TObject.h"
#include "TClonesArray.h"

#include "GcrSelectedXtal.h"
#include "GcrSelectVals.h"

/** @class GcrSelect 
 * @brief Primary Root object containing GcrSelect data.
 *
 *
 * $Header$
*/


/** @page gcrSelectRootData_doc
    @section GcrSelect GcrSelect Documentation
    Notes 

<table>
<tr><th><EM> Method <th> Return Type </EM><th> Description <th>

<small>
<tr><td> getGcrSelectedXtalCol
<td>const TClonesArray*<td> Returns a const pointer to the TClonesArray of 
GcrSelectedXtal
<tr><td> getGcrSelectVals
<td>const GcrSelectVals*<td> Returns const pointer to GcrSelectVals object

</small>
</table>
*/

/** @page gcrSelectRootData_doc
    @section GcrSelect GcrSelect Documentation
 

<table>
<tr>
<th>Method</th><th>Description</th><th>Return Value</th>
</tr>

<small>
<tr>
<td>getGcrSelectedXtalCol( )</td><td>Returns the Collection of hits selected for calibration</td><td>TClonesArray* of SelectedXtals</td>
</tr><tr>
<td>getGcrSelectVals( )</td><td>Returns a set of simple parameters for selected events</td><td>GcrSelectVals*</td>
</tr>

</small>
</table>

*/



class GcrSelect : public TObject
{

public:

    GcrSelect();

    virtual ~GcrSelect();

    void initialize();

    void Clear(Option_t *option="");
    void Fake( Int_t ievent, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const GcrSelect &, const std::string & name = "" ) const; // for tests
    void Print(Option_t *option="") const;
    const TClonesArray* getGcrSelectedXtalCol() const { return m_gcrSelectedXtalCol; };
    const GcrSelectVals* getGcrSelectVals() const { return m_gcrSelectVals; };
    
    /// C.L. 08/22/06:
    //void addGcrSelectedXtal(GcrSelectedXtal* gcrSelectedXtal){ m_gcrSelectedXtalCol->Add(gcrSelectedXtal); };
    GcrSelectedXtal* addGcrSelectedXtal();
    void addGcrSelectVals(GcrSelectVals* gcrSelectVals){ m_gcrSelectVals= gcrSelectVals; };
   
private:
    
    Int_t m_indGcrSelectXtal;
    TClonesArray* m_gcrSelectedXtalCol;
    GcrSelectVals* m_gcrSelectVals;
    

    ClassDef(GcrSelect,6)
};

#endif
