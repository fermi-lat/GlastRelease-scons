set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=GraphicsSvc -version=v1 -path=d:/packages/gaudi_dev $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

