sampleMinuitFit.c

This is a small sample Root script that displays the use of the TMinuit Root class.  This class is used to do Minuit fitting; in this case it fits four points (with differing y error).  To run the macro, start root and do

.L sampleMinuitFit.c   // Loads the macro.
DoFit()                // Does the fi, using Minuit.
DisplayFit()           // Display the fit line and the points.

Thomas Lindner - March 7, 2001

