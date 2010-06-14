set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=reconRoot -version=v0 -path=I:/packages/Root2IDL $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

