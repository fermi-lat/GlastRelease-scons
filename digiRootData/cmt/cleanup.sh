tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=digiRootData -version=v0 -path=I:/packages/GlastGaudi $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

