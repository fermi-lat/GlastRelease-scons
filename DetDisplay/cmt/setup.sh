# echo "Setting DetDisp v1 in c:/packages/gui_dev"

CMTROOT=c:\packages/CMT/v1r5p1; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=DetDisp -version=v1 -path=c:/packages/gui_dev $* >$tempfile; . $tempfile; /bin/rm -f $tempfile
