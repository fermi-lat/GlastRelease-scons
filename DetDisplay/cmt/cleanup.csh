set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=DetDisp -version=v1 -path=c:/packages/gui_dev $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

