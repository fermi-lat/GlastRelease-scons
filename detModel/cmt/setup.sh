# echo "Setting detModel v1r1 in /scratch/users/riccardo/glast-slac"

CMTROOT=/atl2lx/users/riccardo/CMT/v1r6; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=detModel -version=v1r1 -path=/scratch/users/riccardo/glast-slac $* >${tempfile}; . ${tempfile}; /bin/rm -f ${tempfile}
