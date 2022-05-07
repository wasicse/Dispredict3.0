#!/bin/bash

DIRNAME=$1

#change to your own path here
SCRIPT=$(readlink -f "$0")
SHPATH=$(dirname "$SCRIPT")
if [ $# -eq 0 ];then
	echo "Usage: ./DisoComb.sh fastafile"
	exit
fi

if [ ! -e $DIRNAME ];then
	echo "File don't exist\n"
fi

if [ "${DIRNAME:0:1}" = "/" ];then
    ABSDIR=`dirname $DIRNAME`
else
    ABSDIR="`pwd`"/"`dirname $DIRNAME`"
fi

FILENAME=${DIRNAME##*/}
FNAME=${FILENAME%.*}
TMPDIR=$ABSDIR/"tmp_"$FNAME
SCOREDIR=$ABSDIR/"features"_$FNAME
PREDIR=$ABSDIR/"pred"_$FNAME

if [ ! -d $TMPDIR ];then
	mkdir $TMPDIR
fi

if [ ! -d $SCOREDIR ];then
	mkdir $SCOREDIR
fi

if [ ! -d $PREDIR ];then
	mkdir $PREDIR
fi


cd $SHPATH/programs/DisoRDPbind/
./DisoRDPbind $ABSDIR/$FILENAME $TMPDIR/${FNAME}_disordpbind.predictions
cd $SHPATH/programs/fMoRFpred/
./fMoRFpred.sh $ABSDIR/$FILENAME $TMPDIR/${FNAME}_fmorfpred.predictions
cd $SHPATH/programs/DFLpred/
java -jar DFLpred.jar $ABSDIR/$FILENAME $TMPDIR/${FNAME}_dflpred.predictions

ls $TMPDIR/*.seq|while read line;
do 
	a=${line##*/}
	id=${a%%.seq}
	echo $id
done > $TMPDIR/idlist

cd $SHPATH/programs/blast-2.2.24
cat $TMPDIR/idlist|while read id;
do
	./bin/psiblast -query $TMPDIR/$id.seq -db ./db/swissprot -num_iterations 3 -out $TMPDIR/$id.out -out_ascii_pssm $TMPDIR/$id.pssm 
	if [ ! -e  $TMPDIR/$id.pssm ];then
		touch $PREDIR/use_default_pssm_$id
		$SHPATH/programs/create_default_pssm  $TMPDIR/$id.seq > $TMPDIR/$id.pssm
	fi
done

export IUPred_PATH=$SHPATH/programs/iupred
cd $SHPATH/programs/iupred
cat $TMPDIR/idlist|while read id
do
	./iupred $TMPDIR/$id.seq long |grep -v '^#' > $TMPDIR/$id.long
	./iupred $TMPDIR/$id.seq short |grep -v '^#' > $TMPDIR/$id.short
done

#generate PL from prediction by DFLpred
cat $TMPDIR/idlist|while read line;do
	grep '>'$line -A 2 $TMPDIR/${FNAME}_dflpred.predictions|sed 's/,/ /g' > $TMPDIR/$line.dfl
done

#generate PD PR PP from prediction by DisoRDPbind
cat $TMPDIR/idlist|while read line;do
	grep '>'$line -A 7 $TMPDIR/${FNAME}_disordpbind.predictions|sed 's/,/ /g' |sed "s/.*binding.*://g"  > $TMPDIR/$line.rdp
done

#generate PM from prediction by fMoRFpred
cat $TMPDIR/idlist|while read line;do
	grep '>'$line -A 4 $TMPDIR/${FNAME}_fmorfpred.predictions|sed 's/,/ /g' > $TMPDIR/$line.fmorf
done

ls $TMPDIR/*.pssm|while read line;
do
	sed '1,3d' $line|tac|sed '1,6d'|tac  > tmp
	mv $line $line.bak
	mv tmp $line
done 

cd $SHPATH/programs/logReg/
cat $TMPDIR/idlist|while read line;do
	./logitReg $TMPDIR $line 
	mv $TMPDIR/$line.score $SCOREDIR/$line.score
	mv $TMPDIR/$line.log.pred $PREDIR/$line.log.pred
done

cd $SHPATH/programs/NNpackage/
number=$RANDOM
mkdir tmp$number
cat $TMPDIR/idlist|while read line;do
	cut -d $'\t' -f 3-317 $SCOREDIR/$line.score > ./tmp$number/$line.ttscore
	cut -d $'\t' -f 1-2 $SCOREDIR/$line.score > ./tmp$number/$line.ttindex
	python3 Disnet.py ./tmp$number/$line.ttscore > ./tmp$number/$line.ttpreds
	paste ./tmp$number/$line.ttindex ./tmp$number/$line.ttpreds > $PREDIR/$line.nn.pred
	cp -r ./tmp$number/* ../../output/

	if [ -d tmp$number ];then
		rm tmp$number/$line.*
	fi
done

if [ -d tmp$number ];then
	rmdir tmp$number
fi


if [ -d $TMPDIR ];then
	rm $TMPDIR/*
	rmdir $TMPDIR
fi

if [ -d $SCOREDIR ];then
	rm $SCOREDIR/*
	rmdir $SCOREDIR
fi

if [ -d $ABSDIR/$FNAME ];then
	rm $ABSDIR/$FNAME/*
	rmdir $ABSDIR/$FNAME	
fi


