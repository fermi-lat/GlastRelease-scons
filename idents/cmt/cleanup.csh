set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=idents -version=v1 -path=d:/packages/gaudi_dev $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

