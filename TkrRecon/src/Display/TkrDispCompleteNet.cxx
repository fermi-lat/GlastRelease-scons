//------------------------------------------------------------------------------
// TkrDispCompleteNet for TkrNeurlPR Implementation
//
// Used for displaying the entire neural network. 
// Mainly for debugging.  Copied from Tracy's display routines.
//
// b. allgood and w. atwood, 3/02  
//------------------------------------------------------------------------------

#include "TkrRecon/Display/TkrDispCompleteNet.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

TkrDispCompleteNet::TkrDispCompleteNet(IDataProviderSvc* dataProviderSvc, 
                                       ITkrGeometrySvc* pTkrGeometry)
{
    dps     = dataProviderSvc;
    pTkrGeo = pTkrGeometry;
}

//-------------------- private ----------------------

void TkrDispCompleteNet::update()
{
    TkrCandidates* pTkrCandidates = SmartDataPtr<TkrCandidates>(dps,"/Event/TkrRecon/TkrCandidates");

	//Now see if we can do the drawing
	if (pTkrCandidates)
	{

		TkrNeuralNet* pTkrNeuralNet = dynamic_cast<TkrNeuralNet*>(pTkrCandidates);

        int numDispNeurons = pTkrNeuralNet->numNeurons();
        int colorIdx      = 6;
		const char violet[]     = "violet";
		setColor(violet);

        gui::DisplayRep* pDisplay = this;

		TkrNeuralNet::TkrNeuronList tmpList = pTkrNeuralNet->neurons();
		TkrNeuralNet::TkrNeuronList::const_iterator hypo;

		for(hypo  = tmpList.begin(); 
		hypo != tmpList.end();   hypo++){
  
			Point point0 = (*hypo).getPnt(top);
			Point point1 = (*hypo).getPnt(bottom);

			moveTo(point0);
			lineTo(point1);

		}

	    setColor("blue");
    }

    return;
}