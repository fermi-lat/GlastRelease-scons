tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=GraphicsSvc -version=v1 -path=d:/packages/gaudi_dev $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

