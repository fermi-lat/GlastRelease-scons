set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=GuiSvc -version=v1r3 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/pdrApp_v2 $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

