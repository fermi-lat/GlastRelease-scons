tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=DetDisp -version=v1 -path=c:/packages/gui_dev $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

