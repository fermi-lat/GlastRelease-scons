tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=detModel -version=v1 -path=/afs/slac.stanford.edu/u/ey/jrb/glast/jrbPack $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

