#!/bin/bash

dir=~/data/decoys/bqian/t196/
list=tmp.list
log=~/decoyfeatures/test.log
prot=t196
paths=paths.txt.gsbr

whips=(whip01 whip02 whip03 whip04 whip05 whip06 whip07 whip08 whip09 whip10 whip11)
whip=0

function rundecoyfeatures
{
    dir=$1; list=$2; log=$3; prot=$4; paths=$5

    args='-decoyfeatures features -score -fa_input -try_both_his_tautomers -constant_seed'
    cmd="nice ~/rosetta_features.gcc df $prot _ -paths ${paths} ${args} -l ${list}"
    ssh -x ${whips[whip]} "cd $dir; rm -f df${prot}.fasc; $cmd > ${log}; " &
    tmp=`expr $whip + 1`
    whip=`expr $tmp % 11`
}

#rundecoystats ~/data/decoys/bqian/t196/ tmp.list ~/tmp.log t196 paths.txt

function rundf {
    dirh=~/data/decoys/
    logh=~/decoyfeatures/log/
    
    name=bqian
    for i in t196  t199  t205  t206  t223  t232  t234  t249  t262  t279; do
	dir="~/data/decoys/${name}/$i"
	for job in gsbr lowscore lowrms relax_native; do
	    echo $name $i
	    rundecoyfeatures $dir $job.list ${logh}$job.$name.$i.df.log $i paths.txt.$job
	done
    done
    
    name=kira
    for i in 1aa3 1acp 1ail 1b72 1bf4 1erv 1ig5 1pgx 1r69 1tif 1tig 1utg 1vii 256b; do
	dir="~/data/decoys/${name}/$i"
	for job in gsbr lowscore lowrms relax_native; do
	    echo $name $i
	    rundecoyfeatures $dir $job.list ${logh}$job.$name.$i.df.log $i paths.txt.$job
	done
    done
    
    name=sraman
    for i in 1fvk  1hb6  1hoe  1who  3il8  4rnt; do
	dir="~/data/decoys/${name}/$i"
	for job in gsbr lowscore lowrms relax_native; do
	    echo $name $i
	    rundecoyfeatures $dir $job.list ${logh}$job.$name.$i.df.log $i paths.txt.$job
	done
    done
}

#############################################

function listdf 
{
#    tmp=`pwd`
    cd ~/data/decoys/
    for i in `'ls'`; do
	cd ~/data/decoys/
	if [ -d $i ]; then
	    cd ~/data/decoys/$i
	#echo $i `'ls'`
	    for j in `'ls'`; do
		cd ~/data/decoys/$i/$j
		for k in relax_native gsbr lowscore lowrms; do
		    echo $i $j $k `'ls' out/$k | wc -l `
		done
	    #ls out
	    #mkdir out
	    #mkdir out/relax_native
	    #mkdir out/gsbr
	    #mkdir out/lowscore
	    #mkdir out/lowrms
	    #cp ../../paths.txt.*.template .	    
	    #for job in gsbr lowscore lowrms relax_native; do
	#	mv paths.txt.$job.template paths.txt.$job
	#    done
	    done
	fi
#	cd $tmp
    done
}

#############################################

function whipkillall 
{
for i in 01 02 03 04 05 06 07 08 09 10 11
do
  echo $i 
  ssh whip$i killall rosetta_energies.gcc
  ssh whip$i killall rosetta_features.gcc
  #ssh whip$i ps aux | grep sheffler | grep rosetta
done
}

function whiplistproc 
{
    tot=0
    for i in 01 02 03 04 05 06 07 08 09 10 11
      do
      tmp=`ssh whip$i 'ps aux | grep rosetta | grep sheffler | grep -v grep | wc -l'`
      echo $i $tmp
      tot=`expr $tot + $tmp`
    done
    echo $to
}

cd ~/decoy_buried_polar/
for i in `/bin/ls log`; do
    grep DF_POLAR log/$i > df_polar/$i.dfpolar    
done

rm -f relax_native_buried_polar.data decoy_bruied_polar.data

for i in `/bin/ls df_polar/*native*`; do
    cat $i >> relax_native_buried_polar.data
done

for i in `/bin/ls df_polar/*lowscore*`; do
    cat $i >> decoy_bruied_polar.data
done

for i in `/bin/ls df_polar/*gsbr*`; do
    cat $i >> decoy_bruied_polar.data
done

###################3
# pdbs
#  -decoyfeatures -score -fa_input -try_both_his_tautomers -chain_inc -constant_seed -skip_missing_residues -l pdb.list6
###################

# do on whip
cd ~/pdb_buried_polar
split -d -a1 -l 400 ~/data/pdb/cullpdb_pc50_res2.0_R0.25_d050604_chains3866.list pdb.list


whips=(whip01 whip02 whip03 whip04 whip05 whip06 whip07 whip08 whip09 whip11)
whip=0
function runpdb
{
    dir=$1; list=$2; log=$3; job=$4
    args="-decoyfeatures $job -score -fa_input -try_both_his_tautomers  -constant_seed"
    cmd="nice ~/rosetta_features.gcc ${args} -l ${list}"
    ssh -x ${whips[whip]} "cd $dir; $cmd > ${log}; " &
    tmp=`expr $whip + 1`
    whip=`expr $tmp % 11`
}

for i in 0 1 2 3 4 5 6 7 8 9; do
    runpdb ~/pdb_buried_polar/ pdb.list$i pdb_buried_polar_with_ss_$i.log polar
done

for i in 0 1 2 3 4 5 6 7 8 9; do
    cat pdb_energies_$i.log >> pdb_buried_polar.log
done

#######################################################3


function copydf {
    name=$1
    prot=$2
    type=$3
    source="/users/sheffler/data/decoys/$name/$prot/out/$type"
    target="/users/sheffler/decoyfeatures/$name/$prot/$type" 
    echo $source $target
    mkdir -p $target
    cp -r $source $target
}

name=bqian
for prot in t196  t199  t205  t206  t223  t232  t234  t249  t262  t279; do
    for type in gsbr lowscore lowrms relax_native; do
	copydf $name $prot $type
    done
done

name=kira
for prot in 1aa3 1acp 1ail 1b72 1bf4 1erv 1ig5 1pgx 1r69 1tif 1tig 1utg 1vii 256b; do
    for type in gsbr lowscore lowrms relax_native; do
	copydf $name $prot $type	
    done
done

name=sraman
for prot in 1fvk  1hb6  1hoe  1who  3il8  4rnt; do
    for type in gsbr lowscore lowrms relax_native; do
	copydf $name $prot $type
    done
done
