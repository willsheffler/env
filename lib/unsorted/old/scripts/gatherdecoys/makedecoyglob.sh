
function makedecoyglob {
    source=$1
    globdirs=$2
    for i in `/bin/ls -F $source | grep /`; do
	if [ -e $source/$i/gather.log ]; then
	    for j in  `/bin/ls -F $source/$i | grep /`; do
		for k in $globdirs; do
		    if [ -d $source/$i/$j/$k ]; then
			#echo $source/$i/$j/$k #`ls $i/$j/$k | wc -l`
			for l in `/bin/ls $source/$i/$j/$k/*.pdb`; do
			    echo $l
			done
		    fi
		done
		#echo `du -sh $source/$i/$j/.RData`
	    done
	fi
    done
}

listfile="$1"
listerr="$1.err"
if [ -z $listfile ]; then
    listfile=/dev/stdout
    listerr=/dev/stderr
fi

source="/users/sheffler/data/decoys/"
globdirs="lowscore lowrms gsgr gsbr rms0 rms1 rms2 rms3 rms4 rms5 rms6 rms7 rms8 rms9"
makedecoyglob $source "$globdirs" > $listfile 2> $listerr
