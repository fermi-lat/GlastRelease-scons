set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=DOGlastHitsEvt -version=v0 -path=D:/code/packages/glast >$tempfile; source $tempfile; /bin/rm -f $tempfile

