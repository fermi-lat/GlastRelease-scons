These macros create several kinds of formats to plot histograms 
from raw and recon data, as follows:

'hit4layers' plots 4 histograms of 4 different HIT Layers, selected by 
the user, in one pad. It outputs a ps file with the result (Raw data).

'hitlayers' loops over ALL Hit Layers pltting in one pad (Raw data). 
It outputs a ps file with the result (Raw data).

'hitstrips' loops over all layers plotting hits per strip histogrmas. 
It creates a ps file with the result (Raw data).

'totmap' loops over all layers plotting tot and creating a ps file with
the result (Raw data).

'epsfilecreator' plots histograms one per pad and creates eps files for 
every of them as selected or commented out by the user (Recon data)
The sequence can be followed to create more histos as the user needs 
them

'reconhistograms' plots all recon histograms created in RootTreeAnalysis
and creates a single ps file in several pages



