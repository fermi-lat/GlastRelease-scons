# echo "Setting detModel HEAD in /scratch/users/riccardo/glast-slac"

CMTROOT=/atl2lx/users/riccardo/CMT/v1r6; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=detModel -version=HEAD -path=/scratch/users/riccardo/glast-slac $* >${tempfile}; . ${tempfile}; /bin/rm -f ${tempfile}
