set tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=GlastEvent -version=v7r2 -path=/scratch/users/frailis/GlastPack $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

