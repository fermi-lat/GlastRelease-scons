tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=GlastEvent -version=v7r2 -path=/scratch/users/frailis/GlastPack $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

