#!/bin/bash

function lsd
{
    for i in `'ls' $1`; do
	if [ -d $1/$i ]; then
	    echo $i
	fi 
    done
}

function gatherdecoys 
{
    srcdir=$1
    target=$2
    natdir=$3
    prot=$4
    dtake=$5
    dleave=$6
    stake=$7
    sleave=$8
    rosdata=$9

    ddirs=`/bin/ls -F $srcdir | grep '/' | egrep $dtake | egrep -v $dleave`
    #echo ddirs $ddirs
    sfile=`/bin/ls $srcdir/*/*.fasc | egrep $stake | egrep -v $sleave`
    #echo sfil $sfile

    if [ -z "$ddirs" ] ; then
	echo "$prot: no decoy dirs which match '$dtake' but not '$dleave'"; 
	return 
    fi
    if [ -z "$sfile" ] ; then
	echo "$prot: no fasc files which match '$stake' but not '$sleave'"; 
	return 
    fi

    wd=`pwd`
    mkdir -p $target/log
    cd $target

    rm -f $target/log/decoy_dirs.txt 
    rm -f $target/log/score_files.txt
    for j in $ddirs;  do 
	echo $srcdir/$j >> $target/log/decoy_dirs.txt ; 
    done
    for j in $sfile; do 
	echo $j >> $target/log/score_files.txt ; 
    done

    echo $prot `wc -l log/score_files.txt` `wc -l log/decoy_dirs.txt`

    echo "     making fasc files (via R script) ..."
    R CMD BATCH ~/scripts/gatherdecoys/make.fasc.files.R log/make.fasc.files.Rout

    for i in `lsd . | egrep -v "relax_native|source|log|rosetta_data"`; do
 	echo "     copying decoys to $i"
 	pdbs=`awk '{ print( $1 ) }' $i.fasc | grep -v filename`
 	#echo $pdbs
 	for j in $pdbs; do
	    if [ -s $j ]; then
		#echo $j > /dev/null
		nice cp $j $i/
	    fi
 	done
 	   /bin/ls $i/* > $i.list
    done

    if [ ! -d $target/rms9 ]; then
	echo "ERROR: something went wrong in make.fasc.files.R script!"
	return
    fi

    cp /users/sheffler/data/decoys/paths.txt.template $target/paths.txt

    # copy rosetta fragments, native, etc....
    if [ -z "`/bin/ls $rosdata/*$prot*`" ]; then 
	echo "     WARNING: no rosetta data found for $prot! in"
	echo "     $rosdata"
    else
	echo "     copying fragments from $rosdata"
	mkdir -p $target/rosetta_data/ 
	nice cp $rosdata/*$prot* $target/rosetta_data/
    fi

    # copy over relax natives, if available, and score then if needed
#     if [ -d $natdir ]; then
# 	echo "     copying relax natives from $protfrom $natdir"
# 	mkdir -p $target/relax_native
# 	nice cp $natdir/* $target/relax_native/
# 	if [ -f $natdir.list ]; then 
# 	    cp $natdir.list $target/ ;
# 	else 
# 	    /bin/ls relax_native/*.pdb > relax_native.list
# 	fi
# 	if [ -f $natdir.fasc ]; then 
# 	    cp $natdir.fasc $target/ ; 
# 	elif [ -d rosetta_data ]; then
# 	    echo "WARNING: no fasc available for relax_natives"
# 	    chain=`/bin/ls $target/rosetta_data | egrep "..$prot.03......._v..3"`
# 	    #echo $chain
# 	    if [ -f rosetta_data/$chain ] && [ -f /users/sheffler/rosetta.gcc ]; then
# 		chain=${chain:6:1}
# 		#echo sc $prot $chain
# 		nice /users/sheffler/rosetta.gcc sc $prot $chain -score -scorefile relax_native -fa_input -constant_seed -try_both_his_tautomers -l relax_native.list > log/score_relax_native.log 2> log/score_relax_native.err
# 	    else
# 		echo "WARNING: can't score relax_natives"
# 	    fi
# 	fi
#     else
# 	echo "     WARNING: no relax natives found for $prot! in"
# 	echo "     $natdir"
#     fi

    cd $wd
}
