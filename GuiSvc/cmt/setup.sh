# echo "Setting GraphicsSvc v1 in d:/packages/gaudi_dev"

CMTROOT=d:/packages/CMT/v1r5p1; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=GraphicsSvc -version=v1 -path=d:/packages/gaudi_dev $* >$tempfile; . $tempfile; /bin/rm -f $tempfile
