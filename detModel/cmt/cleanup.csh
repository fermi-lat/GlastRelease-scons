set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=detModel -version=v1 -path=/afs/slac.stanford.edu/u/ey/jrb/glast/jrbPack $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

