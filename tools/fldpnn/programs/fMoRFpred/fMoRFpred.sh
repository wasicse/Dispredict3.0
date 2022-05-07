#!/bin/sh
java -jar fMoRFpred.jar $1;
cd ${1%.*};
ls -trQ | xargs cat > $2;
cd ..;
#rm ${1%.*}/*;
#rmdir ${1%.*};

