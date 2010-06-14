tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=DOGlastHitsEvt -version=v0 -path=D:/code/packages/glast >$tempfile; . $tempfile; /bin/rm -f $tempfile

