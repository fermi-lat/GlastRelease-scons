# echo "Setting reconRoot v0 in I:/packages/Root2IDL"

CMTROOT=I:\packages/CMT/v1r5p1; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=reconRoot -version=v0 -path=I:/packages/Root2IDL $* >$tempfile; . $tempfile; /bin/rm -f $tempfile
