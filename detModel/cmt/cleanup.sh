tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=detModel -version=HEAD -path=/scratch/users/riccardo/glast-slac $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

