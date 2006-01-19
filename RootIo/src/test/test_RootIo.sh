#!/usr/bin/env bash

#=================================================
# compile
#=================================================

cd ../../../../Gleam/*/cmt
cmt bro make || eval 'echo GLEAM BUILD FAILED ; exit' ;
cd ../../../RootIo/*/cmt
make test || eval 'echo VALIDATION BUILD FAILED ; exit' ;


#=================================================
# write job
#=================================================

cd ../../../Gleam/*/cmt
. setup.sh
cd ../../../rootTestData/*/data/vertical_surface_muons
../../../../Gleam/*/$CMTCONFIG/Gleam.exe ./jobOptions.txt


#=================================================
# read job
#=================================================

cd ../../../../RootIo/*/cmt
. setup.sh
cd ../src/test
rm -f jobOptions.log
../../$CMTCONFIG/test_RootIo.exe jobOptions.txt | grep -E '(Deposited Energy|VolumeID|TDS)' > jobOptions.log
diff -q -s jobOptions.log jobOptions.ref || eval 'echo VALIDATION DIFF FAILED ; exit' ;


#=================================================
# end
#=================================================

echo Everything OK

