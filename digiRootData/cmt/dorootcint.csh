# do rootcint for Digi ROOT classes

# Save the current directory, so we can return when finished
echo "cd" `pwd` > ${DIGIROOTDATAROOT}\changeDir.csh

cd ${DIGIROOTDATAROOT}\digiRootData
# Remove old cint files
rm EventCint.cxx
rm EventCint.h

${ROOT_PATH}\bin\rootcint -f EventCint.cxx -c AcdHeader.h AcdId.h AcdTile.h CalHeader.h CalLog.h ESAPID.h DigiEvent.h L1T.h LiveTime.h LogId.h StripId.h SubSystemHeader.h Tagger.h TkrHeader.h TkrLayer.h EventLinkDef.h

# Return to the original directory
source ${DIGIROOTDATAROOT}\changeDir.csh
rm ${DIGIROOTDATAROOT}\changeDir.csh
