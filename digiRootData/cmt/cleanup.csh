set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=digiRootData -version=v0 -path=I:/packages/GlastGaudi $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

