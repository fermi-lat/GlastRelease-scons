tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=reconRoot -version=v0 -path=I:/packages/Root2IDL $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

